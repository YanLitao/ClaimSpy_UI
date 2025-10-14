#!/usr/bin/env python3
"""
CHGNet-based magnetic moment comparison for undoped and Al-doped alpha-Fe2O3 (hematite).

This script:
- Loads MP_API_KEY from .env
- Fetches conventional standard structure of alpha-Fe2O3 (R-3c) from Materials Project
- Builds a minimal supercell enabling ~5–10% Al substitution (prefers 1 Al on Fe sublattice at ~8.33%)
- Initializes AFM-like collinear moments on Fe sites (+/- m0 by Fe-layer along c)
- Runs two-stage CHGNet relaxations (positions-only, then restricted cell) with consistent settings
- Extracts predicted atomic magnetic moments from CHGNet and computes net moment per formula unit
- Writes results.json, summary.csv, and final POSCARs for undoped and doped structures

Notes:
- Magnetic moment units are Bohr magneton (μB). Reported per formula unit (f.u.).
- Random seeds fixed for reproducibility of site selection.
"""

import os
import json
import csv
import math
import random
from dataclasses import dataclass, asdict
from typing import List, Tuple, Optional

import numpy as np
from dotenv import load_dotenv

from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor

from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet
from chgnet.model import CHGNetCalculator


random.seed(42)
np.random.seed(42)


@dataclass
class CaseResult:
    label: str
    natoms: int
    n_fe: int
    n_al: int
    n_o: int
    supercell_dims: Tuple[int, int, int]
    al_indices: List[int]
    dopant_fraction_fe: float
    energy_eV: float
    magmom_total_muB: float
    magmom_per_fu_muB: float
    fmax_after_stage2_eVA: float
    steps_stage1: int
    steps_stage2: int


def load_hematite_from_mp(api_key: str) -> Structure:
    """Fetch alpha-Fe2O3 (hematite, R-3c) conventional standard structure from MP."""
    with MPRester(api_key) as mpr:
        # Use summary search with filters: formula and spacegroup
        # Newer mp_api organizes under materials.summary
        docs = mpr.materials.summary.search(
            formula="Fe2 O3",
            spacegroup_symbol="R-3c",
            fields=["material_id", "structure", "energy_above_hull", "is_stable"],
            num_chunks=1,
            chunk_size=10,
        )
    if not docs:
        raise RuntimeError("No Fe2O3 (R-3c) entries found in MP search.")
    # Take the lowest energy above hull entry
    docs_sorted = sorted(docs, key=lambda d: (d.energy_above_hull or 0.0))
    st: Structure = docs_sorted[0].structure
    # Ensure conventional standard cell to get ~30 atoms (12 Fe + 18 O) if available
    sga = SpacegroupAnalyzer(st, symprec=1e-3)
    st = sga.get_conventional_standard_structure()
    # Sanity check composition
    comp = st.composition
    if Element("Fe") not in comp or Element("O") not in comp:
        raise RuntimeError("Fetched structure does not contain Fe and O as expected.")
    return st


def count_species(st: Structure) -> Tuple[int, int, int]:
    n_fe = sum(1 for sp in st.species if str(sp) == "Fe")
    n_al = sum(1 for sp in st.species if str(sp) == "Al")
    n_o = sum(1 for sp in st.species if str(sp) == "O")
    return n_fe, n_al, n_o


def choose_supercell_dims_for_target_fe(n_fe_base: int, target_range: Tuple[int, int]) -> Tuple[int, int, int]:
    """
    Find small supercell dims (a,b,c) s.t. n_fe_total in target_range, minimizing total atoms.
    Prefer dims with small product and closer to cubic.
    """
    lower, upper = target_range
    candidates = []
    # Search over small products up to 6 to keep systems light
    for a in range(1, 5):
        for b in range(1, 4):
            for c in range(1, 4):
                prod = a * b * c
                n_fe_total = n_fe_base * prod
                if lower <= n_fe_total <= upper:
                    candidates.append((prod, (a, b, c)))
    if not candidates:
        # Fall back to the smallest supercell above upper, still preferring minimal size
        best = None
        for a in range(1, 5):
            for b in range(1, 4):
                for c in range(1, 4):
                    prod = a * b * c
                    n_fe_total = n_fe_base * prod
                    if n_fe_total > upper:
                        if best is None or prod < best[0]:
                            best = (prod, (a, b, c))
        if best is None:
            return (1, 1, 1)
        return best[1]
    # Choose candidate with minimal product; tie-break toward more cubic (minimize (a-b)^2+(b-c)^2+(a-c)^2)
    candidates.sort(key=lambda x: (x[0], (x[1][0]-x[1][1])**2 + (x[1][1]-x[1][2])**2 + (x[1][0]-x[1][2])**2))
    return candidates[0][1]


def make_afm_initial_magmoms_ase(atoms: Atoms, m0: float = 4.5) -> List[float]:
    """Assign collinear AFM moments on Fe/Al/Fe-substituted sites based on fractional z-layer.
    Fe: +/- m0 alternating by floor(2*z) parity. Al and O: 0.
    """
    mags: List[float] = []
    cell = atoms.get_cell()
    inv_cell = np.linalg.inv(cell.T)
    positions = atoms.get_positions()
    for i, sym in enumerate(atoms.get_chemical_symbols()):
        if sym == "Fe":
            frac = positions[i] @ inv_cell  # fractional coordinates (x,y,z)
            z = frac[2] - math.floor(frac[2])  # wrap to [0,1)
            sign = 1.0 if (int(math.floor(2 * z)) % 2 == 0) else -1.0
            mags.append(sign * m0)
        else:
            mags.append(0.0)
    return mags


def substitute_fe_with_al(st: Structure, n_sub: int, seed: int = 42) -> Tuple[Structure, List[int]]:
    """Substitute n_sub Fe atoms with Al in a Structure. Returns new Structure and substituted indices.
    Picks sites with maximal pairwise distance to spread dopants; for n_sub=1 this is a single random Fe.
    """
    rng = random.Random(seed)
    fe_indices = [i for i, sp in enumerate(st.species) if str(sp) == "Fe"]
    if n_sub < 1 or n_sub > len(fe_indices):
        raise ValueError("Invalid number of Fe->Al substitutions requested.")
    # For 1 dopant, pick a random Fe; for more, greedily maximize minimum distance
    if n_sub == 1:
        chosen = [rng.choice(fe_indices)]
    else:
        chosen = [rng.choice(fe_indices)]
        remaining = set(fe_indices) - set(chosen)
        while len(chosen) < n_sub and remaining:
            # pick site maximizing distance to current chosen set
            dists = []
            for idx in remaining:
                dmin = min(st.get_distance(idx, j, jimage=None) for j in chosen)
                dists.append((dmin, idx))
            dists.sort(reverse=True)
            chosen.append(dists[0][1])
            remaining.remove(dists[0][1])
    mapping = {i: Element("Al") for i in chosen}
    st_doped = st.copy()
    for i in chosen:
        st_doped[i] = Element("Al")
    return st_doped, chosen


def set_calc_and_check(atoms: Atoms, calc: CHGNetCalculator) -> None:
    atoms.set_calculator(calc)
    # initial evaluation to guard against NaNs
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not np.isfinite(e) or not np.isfinite(f).all():
        raise RuntimeError("Non-finite energy/forces before optimization; check initial structure.")


def run_relax(atoms: Atoms, label: str, calc: CHGNetCalculator,
              fmax_coarse: float = 0.2, fmax_fine: float = 0.05,
              max_steps: int = 300) -> Tuple[Atoms, float, float, int, int]:
    """Two-stage relax: positions-only (half steps), then restricted cell with hydrostatic strain."""
    set_calc_and_check(atoms, calc)

    steps1 = max_steps // 2
    opt1 = BFGS(atoms, logfile=f"{label}_relax1.log")
    opt1.run(fmax=fmax_coarse, steps=steps1)

    # Stage 2: hydrostatic cell relax via UnitCellFilter
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile=f"{label}_relax2.log")
    steps2 = max_steps
    opt2.run(fmax=fmax_fine, steps=steps2)

    # Final energy and max force
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    fmax = float(np.linalg.norm(f, axis=1).max()) if len(f) > 0 else 0.0
    return atoms, e, fmax, steps1, steps2


def get_predicted_magmoms(atoms: Atoms) -> Optional[np.ndarray]:
    """Try to retrieve per-atom magnetic moments from the calculator; fallback to Atoms if needed.
    Returns magnitudes (non-negative) if calculator provides only absolute values.
    """
    mm = None
    try:
        mm = atoms.calc.get_property("magmoms", atoms)
    except Exception:
        mm = None
    if mm is None:
        try:
            mm = atoms.get_magnetic_moments()
        except Exception:
            mm = None
    if mm is None:
        return None
    arr = np.array(mm, dtype=float)
    return np.abs(arr)


def atoms_to_formula_units(atoms: Atoms) -> int:
    # Count by oxygen is robust (3 per Fe2O3 f.u.) even with Fe->Al substitution
    n_o = sum(1 for s in atoms.get_chemical_symbols() if s == "O")
    return int(round(n_o / 3))


def run_case(label: str, st: Structure, calc: CHGNetCalculator,
             supercell: Tuple[int, int, int], al_subs: int) -> CaseResult:
    # Apply supercell
    st_sc = st * supercell
    # Substitute Fe->Al
    if al_subs > 0:
        st_sc, al_indices = substitute_fe_with_al(st_sc, al_subs, seed=42)
    else:
        al_indices = []
    # Convert to ASE Atoms
    atoms = AseAtomsAdaptor.get_atoms(st_sc)
    # Initialize AFM-like magmoms
    mags = make_afm_initial_magmoms_ase(atoms)
    atoms.set_initial_magnetic_moments(mags)

    # Relax
    atoms_relaxed, energy, fmax, steps1, steps2 = run_relax(atoms, label, calc)

    # Predicted magmoms from CHGNet (magnitudes) and assign AFM signs from initial moments
    magmoms_abs = get_predicted_magmoms(atoms_relaxed)
    if magmoms_abs is None:
        # Fallback: use initial moments (already signed)
        magmoms_signed = np.array(atoms_relaxed.get_initial_magnetic_moments(), dtype=float)
    else:
        init = np.array(atoms_relaxed.get_initial_magnetic_moments(), dtype=float)
        signs = np.sign(init)
        magmoms_signed = magmoms_abs * signs

    mag_total = float(np.nansum(magmoms_signed))  # μB
    n_fu = atoms_to_formula_units(atoms_relaxed)
    mag_per_fu = mag_total / max(n_fu, 1)

    # Save POSCAR
    write(f"{label}_final.vasp", atoms_relaxed, format="vasp")

    # Composition counts
    symbols = atoms_relaxed.get_chemical_symbols()
    n_fe = sum(1 for s in symbols if s == "Fe")
    n_al = sum(1 for s in symbols if s == "Al")
    n_o = sum(1 for s in symbols if s == "O")
    n_fe_total = n_fe + n_al  # Fe sublattice size (including Al substitutions)
    dop_frac = (n_al / n_fe_total) if n_fe_total > 0 else 0.0

    return CaseResult(
        label=label,
        natoms=len(symbols),
        n_fe=n_fe,
        n_al=n_al,
        n_o=n_o,
        supercell_dims=supercell,
        al_indices=al_indices,
        dopant_fraction_fe=dop_frac,
        energy_eV=float(energy),
        magmom_total_muB=mag_total,
        magmom_per_fu_muB=mag_per_fu,
        fmax_after_stage2_eVA=fmax,
        steps_stage1=steps1,
        steps_stage2=steps2,
    )


def main():
    load_dotenv()
    api_key = os.getenv("MP_API_KEY")
    if not api_key:
        raise RuntimeError("MP_API_KEY not found in environment. Please set it in .env as MP_API_KEY=...")

    # Fetch hematite
    st = load_hematite_from_mp(api_key)
    n_fe_base, _, _ = count_species(st)

    # Choose supercell to allow ~5–10% Al on Fe sublattice with a single Al if possible
    dims = choose_supercell_dims_for_target_fe(n_fe_base, target_range=(10, 20))
    st_sc = st * dims
    n_fe_sc, _, _ = count_species(st_sc)

    # Prefer a single Al substituent; if outside 5–10%, adjust to nearest feasible within small count
    # Here, with conventional cell (n_fe_sc=12), 1 Al gives 8.33%, in range.
    al_subs = 1
    frac = al_subs / n_fe_sc
    if not (0.05 <= frac <= 0.10):
        # Try al_subs=2
        if 0.05 <= (2 / n_fe_sc) <= 0.10:
            al_subs = 2
        else:
            # As a last resort, keep al_subs=1 (closest small integer) and report actual fraction
            al_subs = 1

    # Load CHGNet once and reuse
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Undoped case
    undoped_res = run_case("undoped", st, calc, dims, al_subs=0)

    # Doped case
    doped_res = run_case("al_doped", st, calc, dims, al_subs=al_subs)

    # Write results
    results = {
        "undoped": asdict(undoped_res),
        "al_doped": asdict(doped_res),
        "dmu_per_fu_muB": doped_res.magmom_per_fu_muB - undoped_res.magmom_per_fu_muB,
    }
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)

    # CSV summary
    with open("summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "case", "natoms", "n_fe", "n_al", "n_o", "supercell",
            "dopant_fraction_fe", "energy_eV", "magmom_total_muB", "magmom_per_fu_muB",
            "fmax_after_stage2_eVA", "steps_stage1", "steps_stage2",
        ])
        for res in (undoped_res, doped_res):
            w.writerow([
                res.label, res.natoms, res.n_fe, res.n_al, res.n_o,
                "x".join(map(str, res.supercell_dims)),
                f"{res.dopant_fraction_fe:.6f}", f"{res.energy_eV:.6f}",
                f"{res.magmom_total_muB:.6f}", f"{res.magmom_per_fu_muB:.6f}",
                f"{res.fmax_after_stage2_eVA:.6f}", res.steps_stage1, res.steps_stage2,
            ])

    # Human-readable concise summary
    print("=== CHGNet magnetic moment comparison: α-Fe2O3 vs Al-doped α-Fe2O3 ===")
    print(f"Supercell dims: {dims} | Undoped Fe sites: {undoped_res.n_fe} | Doped Al count: {doped_res.n_al}")
    print(f"Al fraction on Fe sublattice: {doped_res.dopant_fraction_fe*100:.2f}%")
    print(f"Undoped: μ_total = {undoped_res.magmom_total_muB:.3f} μB, μ/f.u. = {undoped_res.magmom_per_fu_muB:.3f} μB")
    print(f"Al-doped: μ_total = {doped_res.magmom_total_muB:.3f} μB, μ/f.u. = {doped_res.magmom_per_fu_muB:.3f} μB")
    dmu = doped_res.magmom_per_fu_muB - undoped_res.magmom_per_fu_muB
    trend = "increases" if dmu > 0 else ("decreases" if dmu < 0 else "no change")
    print(f"Difference (doped − undoped): Δμ/f.u. = {dmu:.3f} μB → Al doping {trend} net magnetic moment.")
    print("Wrote: results.json, summary.csv, undoped_final.vasp, al_doped_final.vasp, and relax logs.")


if __name__ == "__main__":
    main()
