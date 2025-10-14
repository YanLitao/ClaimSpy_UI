#!/usr/bin/env python3
"""
CHGNet-based estimation of saturation magnetization for undoped and Co-doped hematite (Fe2O3).

This script:
  1) Loads MP_API_KEY from .env and fetches the ground-state hematite structure from Materials Project.
  2) Builds undoped, ~5% Co and ~10% Co substitutional models (Co on Fe sites) via appropriate supercells.
  3) Performs two-stage CHGNet relaxations (positions-only, then hydrostatic cell relax) with conservative settings.
  4) Predicts site magnetic moments via the CHGNet model and converts to Ms in μB/f.u. and emu/g.
  5) Writes standardized outputs (results.json, summary.csv, final POSCARs) and prints a concise summary.

Assumptions:
  - Saturation magnetization is approximated as the sum of per-site magnetic moment magnitudes (μB) per formula unit
    (i.e., spins aligned). If vector moments are available, their magnitudes are used.
  - For undoped Fe2O3 (hematite), the primitive f.u. is Fe2O3 (5 atoms). Doped structures preserve 5 atoms per f.u.
  - Doping levels are approximate: 1 Co in 60-atom (24 Fe) cell ≈ 4.17% (reported as ~5%); 2 Co in the same cell ≈ 8.33% (~10%).

Execution notes:
  - Uses a single CHGNet model instance and a single CHGNetCalculator instance to avoid repeated loads.
  - Avoids interactive blocking calls; all plots/logs are written to files in CWD.
"""

import os
import json
import csv
import math
import random
from dataclasses import dataclass
from typing import List, Tuple, Dict, Any

import numpy as np

from dotenv import load_dotenv

# Pymatgen / ASE / CHGNet imports
from pymatgen.core import Structure, Element
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator

# Materials Project API (mp-api)
try:
    from mp_api.client import MPRester  # newer mp-api
except Exception:  # pragma: no cover
    from pymatgen.ext.matproj import MPRester  # fallback to legacy if needed


AVOGADRO = 6.02214076e23  # 1/mol
MU_B_EMU = 9.274009994e-21  # 1 Bohr magneton in emu (1 A·m^2 = 1000 emu)
AMU_TO_G = 1.66053906660e-24  # grams per atomic mass unit (approx g per atom where mass is in g/mol)

np.random.seed(42)
random.seed(42)


@dataclass
class CalcResult:
    label: str
    n_atoms: int
    n_fu: float
    composition: Dict[str, float]
    co_fraction_fe_sites: float
    muB_total_abs: float
    muB_per_fu: float
    emu_per_g: float
    energy_eV: float
    fmax_after: float
    steps: int
    magmom_source: str


def load_mp_structure(api_key: str) -> Structure:
    """Fetch ground-state hematite (Fe2O3) from Materials Project and return a conventional cell Structure.
    Preference: R-3c (hematite) with minimal energy above hull.
    """
    assert api_key, "MP_API_KEY is required in the environment (.env)."
    with MPRester(api_key) as mpr:
        # Query summaries to pick the lowest EAH Fe2O3 with R-3c if possible
        fields = [
            "material_id",
            "structure",
            "energy_above_hull",
            "is_stable",
            "symmetry",
            "formula_pretty",
        ]
        # mp-api uses 'formula' (not chemical_formula) on the SummaryRester
        try:
            docs = mpr.summary.search(formula="Fe2O3", fields=fields, num_chunks=1, chunk_size=200)
        except Exception:
            # Older clients may expose a namespaced path
            docs = mpr.materials.summary.search(formula="Fe2O3", fields=fields, num_chunks=1, chunk_size=200)
    if not docs:
        raise RuntimeError("No Fe2O3 documents returned from Materials Project API.")

    # Filter R-3c if present, then choose minimal EAH
    def is_r3c(d):
        try:
            sym = d.symmetry
            if sym and getattr(sym, "symbol", None):
                return sym.symbol.upper().replace(" ", "") in {"R-3C", "R-3C:H", "R-3C:R"}
        except Exception:
            pass
        return False

    r3c_docs = [d for d in docs if is_r3c(d)]
    pool = r3c_docs if r3c_docs else docs
    chosen = sorted(pool, key=lambda d: (d.energy_above_hull if d.energy_above_hull is not None else 1e9))[0]
    structure = chosen.structure  # Pymatgen Structure
    # Ensure we have a conventional-ish cell with reasonable size (hematite often ~30 atoms conventional)
    try:
        structure = structure.get_conventional_standard_structure()
    except Exception:
        pass
    return structure


def build_supercell(structure: Structure, mult: Tuple[int, int, int]) -> Structure:
    mat = structure.copy()
    mat.make_supercell(mult)
    return mat


def get_fe_indices(structure: Structure) -> List[int]:
    return [i for i, sp in enumerate(structure.species) if str(sp) == "Fe"]


def max_separated_pair(structure: Structure, idx_list: List[int]) -> Tuple[int, int]:
    # Choose a pair of indices with maximal Cartesian distance
    best_pair = (idx_list[0], idx_list[-1])
    best_dist = -1.0
    for i in range(len(idx_list)):
        for j in range(i + 1, len(idx_list)):
            d = structure.get_distance(idx_list[i], idx_list[j])
            if d > best_dist:
                best_dist = d
                best_pair = (idx_list[i], idx_list[j])
    return best_pair


def dope_fe_with_co(structure: Structure, n_co: int) -> Tuple[Structure, float]:
    """Substitute n_co Fe sites with Co. Returns (new_structure, co_fraction_on_fe_sites).
    Uses a deterministic but reasonable spacing strategy.
    """
    s = structure.copy()
    fe_idx = get_fe_indices(s)
    n_fe = len(fe_idx)
    if n_co < 0 or n_co > n_fe:
        raise ValueError("n_co must be between 0 and number of Fe sites")
    if n_co == 0:
        return s, 0.0

    # Choose sites: for 1 Co pick the median Fe index; for 2 Co pick maximally separated pair; else spread roughly evenly
    chosen: List[int] = []
    if n_co == 1:
        chosen = [fe_idx[len(fe_idx) // 2]]
    elif n_co == 2:
        i, j = max_separated_pair(s, fe_idx)
        chosen = [i, j]
    else:
        # Evenly spaced selection over the Fe indices for small n
        step = max(1, len(fe_idx) // n_co)
        chosen = [fe_idx[k] for k in range(0, len(fe_idx), step)][:n_co]

    # Apply substitution
    for ci in sorted(chosen, reverse=True):  # reverse to keep indices valid
        s[ci] = "Co"

    co_frac = float(n_co) / float(n_fe) if n_fe > 0 else 0.0
    return s, co_frac


def atoms_min_distance(atoms) -> float:
    from ase.geometry import get_distances
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()
    dmin = np.inf
    for i in range(len(atoms)):
        dists, _ = get_distances(pos[i], pos, cell=cell, pbc=pbc)
        # Ignore self at index i (distance 0)
        dists = dists[dists > 1e-6]
        if len(dists):
            dmin = min(dmin, float(np.min(dists)))
    return float(dmin)


def relax_with_chgnet(atoms, calc: CHGNetCalculator, label: str,
                      fmax_coarse: float = 0.2, fmax_fine: float = 0.05,
                      max_steps: int = 300) -> Tuple[Any, float, int, float]:
    """Two-stage relaxation: positions-only then hydrostatic cell relax. Returns:
    (atoms_final, energy_eV, fmax_after, steps_total)
    """
    atoms.set_calculator(calc)

    # Clip bad geometries by uniform expansion if any pair is too close
    dmin = atoms_min_distance(atoms)
    if not np.isfinite(dmin) or dmin < 0.7:
        scale = 0.7 / max(dmin, 1e-3)
        atoms.set_cell(atoms.get_cell() * (1.05 * scale), scale_atoms=True)

    # Stage 1: positions only
    opt1 = BFGS(atoms, logfile=f"relax1_{label}.log", maxstep=0.2)
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)

    # Stage 2: restricted cell (hydrostatic) + positions
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile=f"relax2_{label}.log", maxstep=0.2)
    opt2.run(fmax=fmax_fine, steps=max_steps)

    # Final checks
    try:
        e = float(atoms.get_potential_energy())
        f = atoms.get_forces()
        fmax = float(np.max(np.linalg.norm(f, axis=1))) if len(f) else float("nan")
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(f"Failed to evaluate final energy/forces for {label}: {exc}")

    steps_total = int(opt1.get_number_of_steps() + opt2.get_number_of_steps())
    return atoms, e, fmax, steps_total


def predict_site_magmoms(chgnet_model: CHGNet, structure: Structure) -> np.ndarray:
    """Attempt to extract per-site magnetic moments (μB) from CHGNet.
    Tries several API entry points for robustness.
    Returns array of shape (N,) or (N,3). If vectors, each row is a 3D moment.
    """
    # Strategy 1: direct structure prediction (preferred; lets CHGNet handle graph building)
    try:
        pred = chgnet_model.predict_structure(structure)  # type: ignore
        for key in ("magmom", "magmoms", "magnetic_moments"):
            if key in pred and pred[key] is not None:
                arr = np.array(pred[key], dtype=float)
                if arr.ndim == 1:
                    return arr
                if arr.ndim == 2:
                    return arr
    except Exception:
        pass

    # Strategy 2: CrystalGraph -> predict_graph (tune cutoff to avoid isolated atoms)
    try:
        from chgnet.graph import CrystalGraph  # type: ignore
        try:
            cg = CrystalGraph(atom_graph_cutoff=8.0, bond_graph_cutoff=6.0)
            g = cg.convert(structure)
        except Exception:
            # Fallback constructors
            cg = CrystalGraph()
            try:
                g = cg.convert(structure)
            except Exception:
                g = CrystalGraph.from_structure(structure)
        pred = chgnet_model.predict_graph(g)
        for key in ("magmom", "magmoms", "magnetic_moments"):
            if key in pred and pred[key] is not None:
                arr = np.array(pred[key], dtype=float)
                if arr.ndim == 1:
                    return arr
                if arr.ndim == 2:
                    return arr
    except Exception:
        pass

    # Strategy 3: graph-based predict_structure via StructureGraph
    try:
        from chgnet.graph import StructureGraph  # type: ignore
        sg = StructureGraph.from_structure(structure)
        pred = chgnet_model.predict_structure(sg)
        for key in ("magmom", "magmoms", "magnetic_moments"):
            if key in pred and pred[key] is not None:
                arr = np.array(pred[key], dtype=float)
                if arr.ndim == 1:
                    return arr
                if arr.ndim == 2:
                    return arr
    except Exception:
        pass

    raise RuntimeError("Could not obtain site magnetic moments from CHGNet; API surface may differ.")


def calc_ms_metrics(magmoms: np.ndarray, structure: Structure) -> Tuple[float, float, float]:
    """Compute Ms figures of merit:
    - muB_total_abs: sum of per-site |μ| in μB for the full structure (saturation approximation)
    - muB_per_fu: μB per Fe2O3 formula unit
    - emu_per_g: emu per gram (for the bulk composition of the structure)
    """
    # If vector moments, use magnitudes
    if magmoms.ndim == 2 and magmoms.shape[1] == 3:
        mags = np.linalg.norm(magmoms, axis=1)
    else:
        mags = np.array(magmoms).reshape(-1)

    muB_total_abs = float(np.sum(np.abs(mags)))

    # Number of Fe2O3 formula units (5 atoms per f.u.)
    n_fu = structure.num_sites / 5.0
    muB_per_fu = muB_total_abs / n_fu

    # Mass of this structure in grams: sum( atomic_mass(g/mol) / N_A )
    comp = structure.composition
    mass_g = 0.0
    for el, amt in comp.items():
        mass_g += float(Element(str(el)).atomic_mass) * amt / AVOGADRO

    emu_total = muB_total_abs * MU_B_EMU
    emu_per_g = emu_total / mass_g

    return muB_total_abs, muB_per_fu, emu_per_g


def spin_only_magmoms(structure: Structure) -> np.ndarray:
    """Fallback: estimate per-site magnetic moments (μB) from nominal high-spin states.
    Fe3+ -> 5 μB; Co3+ -> 4 μB; O2- -> 0 μB.
    This provides a saturation upper bound and preserves relative trends with Co-for-Fe substitution.
    """
    mags: List[float] = []
    for sp in structure.species:
        s = str(sp)
        if s == "Fe":
            mags.append(5.0)
        elif s == "Co":
            mags.append(4.0)
        else:
            mags.append(0.0)
    return np.array(mags, dtype=float)


def write_poscar(structure: Structure, path: str) -> None:
    structure.to(fmt="poscar", filename=path)


def main():
    load_dotenv()
    api_key = os.environ.get("MP_API_KEY", "").strip()
    if not api_key:
        raise SystemExit("MP_API_KEY missing. Please set it in .env as MP_API_KEY=...")

    # 1) Load hematite and prepare structures
    base = load_mp_structure(api_key)
    # Ensure oxygen/iron order remains canonical for readability
    base.sort()

    # Use the fetched cell as undoped reference
    undoped = base.copy()

    # Determine Fe count in base cell to pick a supercell that has 24 Fe sites (60-atom equivalent)
    n_fe_base = len([1 for sp in base.species if str(sp) == "Fe"])
    if n_fe_base <= 0:
        raise SystemExit("Fetched structure has no Fe sites; unexpected for Fe2O3.")
    # Minimal replication factor m such that n_fe_base * m >= 24
    m = max(1, math.ceil(24 / n_fe_base))
    # Choose a reasonable 3D factorization for m
    mult = (m, 1, 1)
    if m == 2:
        mult = (2, 1, 1)
    elif m == 3:
        mult = (3, 1, 1)
    elif m == 4:
        mult = (2, 2, 1)
    elif m == 5:
        mult = (5, 1, 1)
    elif m == 6:
        mult = (3, 2, 1)

    sc_dope = build_supercell(base, mult)

    # ~5% Co (≈4.17%): 1 Co on 24 Fe sites
    co5_struct, co5_frac = dope_fe_with_co(sc_dope, n_co=1)

    # ~10% Co (≈8.33%): 2 Co on 24 Fe sites
    co10_struct, co10_frac = dope_fe_with_co(sc_dope, n_co=2)

    # 2) CHGNet model and calculator (single instances)
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    adaptor = AseAtomsAdaptor()

    # 3) Relax structures
    results: List[CalcResult] = []
    structures = [
        ("Fe2O3_undoped", undoped),
        ("Fe2O3_Co_x0p0417", co5_struct),
        ("Fe2O3_Co_x0p0833", co10_struct),
    ]

    created_files: List[str] = []

    for label, struct in structures:
        atoms = adaptor.get_atoms(struct)
        atoms, energy_eV, fmax_after, steps = relax_with_chgnet(atoms, calc, label)

        # Update structure from relaxed Atoms
        struct_relaxed = adaptor.get_structure(atoms)

        # Predict site magnetic moments
        magmom_source = "chgnet"
        try:
            magmoms = predict_site_magmoms(chgnet, struct_relaxed)
        except Exception:
            # Fallback to spin-only estimate if CHGNet magnetic moments are unavailable in this environment
            magmoms = spin_only_magmoms(struct_relaxed)
            magmom_source = "spin_only"

        # Compute Ms metrics
        muB_total_abs, muB_per_fu, emu_per_g = calc_ms_metrics(magmoms, struct_relaxed)

        # Co fraction on Fe sublattice
        n_fe_sites = sum(1 for sp in struct_relaxed.species if str(sp) in ("Fe", "Co"))
        n_co_sites = sum(1 for sp in struct_relaxed.species if str(sp) == "Co")
        co_frac = (n_co_sites / n_fe_sites) if n_fe_sites > 0 else 0.0

        # POSCAR
        poscar_name = f"POSCAR_{label}.vasp"
        write_poscar(struct_relaxed, poscar_name)
        created_files.append(poscar_name)

        # Save log files already created by BFGS
        created_files.append(f"relax1_{label}.log")
        created_files.append(f"relax2_{label}.log")

        # Record result
        comp_dict = {el.symbol if isinstance(el, Element) else str(el): float(amt)
                     for el, amt in struct_relaxed.composition.items()}
        n_fu = struct_relaxed.num_sites / 5.0
        results.append(CalcResult(
            label=label,
            n_atoms=struct_relaxed.num_sites,
            n_fu=float(n_fu),
            composition=comp_dict,
            co_fraction_fe_sites=float(co_frac),
            muB_total_abs=float(muB_total_abs),
            muB_per_fu=float(muB_per_fu),
            emu_per_g=float(emu_per_g),
            energy_eV=float(energy_eV),
            fmax_after=float(fmax_after),
            steps=int(steps),
            magmom_source=magmom_source,
        ))

    # 4) Save results
    results_json = "results.json"
    with open(results_json, "w") as f:
        json.dump([r.__dict__ for r in results], f, indent=2)
    created_files.append(results_json)

    summary_csv = "summary.csv"
    with open(summary_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "label", "n_atoms", "n_fu", "Co_frac_Fe_sublattice",
            "muB_total_abs", "muB_per_fu", "emu_per_g",
            "energy_eV", "fmax_after", "steps", "magmom_source"
        ])
        for r in results:
            w.writerow([
                r.label, r.n_atoms, f"{r.n_fu:.2f}", f"{r.co_fraction_fe_sites:.4f}",
                f"{r.muB_total_abs:.4f}", f"{r.muB_per_fu:.4f}", f"{r.emu_per_g:.4e}",
                f"{r.energy_eV:.6f}", f"{r.fmax_after:.4f}", r.steps, r.magmom_source,
            ])
    created_files.append(summary_csv)

    # 5) Human-readable summary
    # Compare relative increases vs undoped
    base_muB_fu = next(r.muB_per_fu for r in results if r.label == "Fe2O3_undoped")
    base_emu_g = next(r.emu_per_g for r in results if r.label == "Fe2O3_undoped")

    print("\nCHGNet Saturation Magnetization Estimates (sums of |μ|):")
    for r in results:
        rel_mu = (r.muB_per_fu / base_muB_fu - 1.0) * 100.0 if base_muB_fu > 0 else float("nan")
        rel_emu = (r.emu_per_g / base_emu_g - 1.0) * 100.0 if base_emu_g > 0 else float("nan")
        print(
            f"  {r.label} | Co on Fe: {r.co_fraction_fe_sites*100:.2f}% | "
            f"μB/f.u.: {r.muB_per_fu:.3f} ({rel_mu:+.1f}% vs undoped) | "
            f"emu/g: {r.emu_per_g:.3e} ({rel_emu:+.1f}% vs undoped)"
        )

    # Write created files manifest for the harness to collect if desired
    manifest_path = "intermediate_files.json"
    with open(manifest_path, "w") as f:
        json.dump(sorted(list(dict.fromkeys(created_files))), f, indent=2)
    print(f"\nWrote outputs: {', '.join(sorted(list(dict.fromkeys(created_files))))}")


if __name__ == "__main__":
    main()
