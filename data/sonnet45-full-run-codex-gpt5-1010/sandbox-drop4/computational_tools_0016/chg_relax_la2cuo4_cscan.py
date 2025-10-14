#!/usr/bin/env python3
"""
CHGNet relaxations for La2CuO4 vs c-axis expansion.

Workflow:
 1) Load CHGNet model once.
 2) Fetch experimental-like tetragonal La2CuO4 (I4/mmm) structure from Materials Project (uses MP_API_KEY from .env). Fallback to an approximate K2NiF4 prototype if unavailable.
 3) Convert to ASE Atoms once.
 4) For c-scale factors [1.00, 1.02, 1.05, 1.10]:
      - Scale only the c lattice vector while keeping a,b fixed.
      - Positions-only relaxation in two stages (BFGS): fmax 0.2 then 0.05 eV/Å.
      - Record: total energy/atom (eV/atom), max force (eV/Å), Cu–O octahedral bond lengths (Å), and whether octahedral geometry is preserved (boolean).
      - Write final POSCAR for each configuration.
 5) Emit results.json and summary.csv and print a concise text summary.

Notes:
 - Cell DOFs are held fixed; only internal coordinates relax.
 - Octahedral preservation criterion: for every Cu, the 6th-nearest O distance ≤ 3.0 Å and a clear 4+2 split exists (d5 - d4 ≥ 0.05 Å).
 - This script is deterministic (random seeds set), loads CHGNet once, and avoids interactive calls.
"""

from __future__ import annotations

import json
import math
import os
import random
import sys
import warnings
from dataclasses import asdict, dataclass
from typing import Dict, List, Tuple

import numpy as np

# Set seeds for determinism
random.seed(42)
np.random.seed(42)

# Load MP API key from .env if available
try:
    from dotenv import load_dotenv  # type: ignore
    load_dotenv()
except Exception:
    pass


def _fetch_la2cuo4_from_mp() -> "Structure | None":
    """Try to fetch a tetragonal La2CuO4 (I4/mmm) structure from Materials Project.

    Returns None if unavailable or if mp_api is not installed.
    """
    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        return None

    # Try the new mp_api client first
    try:
        from mp_api.client import MPRester  # type: ignore

        with MPRester(api_key) as mpr:
            # Request minimal fields to avoid heavy payloads
            docs = mpr.summary.search(
                formula="La2CuO4",
                fields=["material_id", "structure", "spacegroup"],
            )
            # Filter to tetragonal I4/mmm if possible
            candidates = []
            for d in docs:
                try:
                    spg_symbol = getattr(d.spacegroup, "symbol", None)
                except Exception:
                    spg_symbol = None
                if spg_symbol is None:
                    try:
                        spg_symbol = d["spacegroup"]["symbol"]  # type: ignore[index]
                    except Exception:
                        spg_symbol = None
                if spg_symbol and spg_symbol.upper() == "I4/MMM":
                    candidates.append(d)
            if not candidates:
                # If no filtered candidates, fall back to first returned structure
                candidates = list(docs)
            # Sort deterministically by material_id string
            candidates = sorted(
                candidates,
                key=lambda x: str(getattr(x, "material_id", getattr(x, "material_id", "zzz"))),
            )
            for d in candidates:
                s = getattr(d, "structure", None)
                if s is not None:
                    return s
    except Exception:
        pass

    # Try older pymatgen MPRester as a fallback
    try:
        from pymatgen.ext.matproj import MPRester  # type: ignore

        with MPRester(api_key) as mpr:
            entries = mpr.query(
                {"formula_pretty": "La2CuO4"},
                properties=["material_id", "spacegroup", "structure"],
            )
            # Prefer I4/mmm
            pref = [e for e in entries if e.get("spacegroup", {}).get("symbol") == "I4/mmm"]
            if not pref:
                pref = entries
            # Sort by material_id and return first
            pref = sorted(pref, key=lambda e: str(e.get("material_id", "zzz")))
            if pref:
                return pref[0]["structure"]
    except Exception:
        pass

    return None


def _approx_k2nif4_la2cuo4(a: float = 3.80, c: float = 13.20):
    """Construct an approximate I4/mmm La2CuO4 (K2NiF4-type) cell.

    Fractional internal parameters (z) chosen from typical literature values
    for the high-temperature T-phase: z_La ≈ 0.362, z_Oap ≈ 0.175.
    This is only used if MP retrieval is unavailable.
    Z (number of formula units per cell) = 2; total 14 atoms.
    """
    from pymatgen.core import Lattice, Structure

    lat = Lattice.tetragonal(a, c)

    z_La = 0.362
    z_Oap = 0.175

    species = []
    frac_coords = []

    # Cu at 2a: (0,0,0); (0,0,1/2)
    for z in [0.0, 0.5]:
        species.append("Cu")
        frac_coords.append([0, 0, z])

    # La at 4e: (0,0,z_La); (0,0,1-z_La); + (add 1/2 in z for the other pair)
    for z in [z_La, 1 - z_La, z_La + 0.5, 1 - z_La + 0.5]:
        species.append("La")
        frac_coords.append([0, 0, z % 1.0])

    # O planar (O1) at 4c: (0,1/2,0); (1/2,0,0); plus those translated by 1/2 in z
    for z in [0.0, 0.5]:
        species.extend(["O", "O"])
        frac_coords.extend([[0, 0.5, z], [0.5, 0, z]])

    # O apical (O2) at 4e: (0,0,z_Oap); (0,0,1-z_Oap); + (add 1/2 in z)
    for z in [z_Oap, 1 - z_Oap, z_Oap + 0.5, 1 - z_Oap + 0.5]:
        species.append("O")
        frac_coords.append([0, 0, z % 1.0])

    struct = Structure(lat, species, frac_coords, to_unit_cell=True, coords_are_cartesian=False)
    return struct


def get_initial_structure():
    """Return a pymatgen Structure representing tetragonal La2CuO4."""
    s = _fetch_la2cuo4_from_mp()
    if s is not None:
        return s
    # Fallback
    return _approx_k2nif4_la2cuo4()


def atoms_copy_with_c_scale(atoms, c_scale: float):
    """Return a copy of ASE Atoms with only the c lattice vector scaled by c_scale."""
    from ase import Atoms

    a2 = atoms.copy()
    cell = a2.cell.array.copy()
    cell[2, :] *= c_scale
    a2.set_cell(cell, scale_atoms=True)
    return a2


def ensure_finite_energy_forces(atoms, tag: str = ""):
    import numpy as _np

    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not (_np.isfinite(e) and _np.isfinite(f).all()):
        raise RuntimeError(f"Non-finite energy/forces encountered {tag}")


def relax_positions_only(atoms, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300, log_prefix="relax"):
    from ase.optimize import BFGS

    # Stage 1: coarse forces
    opt1 = BFGS(atoms, logfile=f"{log_prefix}_stage1.log")
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)

    # Stage 2: fine forces
    opt2 = BFGS(atoms, logfile=f"{log_prefix}_stage2.log")
    opt2.run(fmax=fmax_fine, steps=max_steps)

    return atoms


def compute_max_force(atoms) -> float:
    f = atoms.get_forces()
    norms = np.linalg.norm(f, axis=1)
    return float(np.max(norms))


def cu_o_bonds(atoms) -> Tuple[List[List[float]], bool]:
    """Compute 6 nearest Cu–O distances for each Cu in Atoms.

    Returns (list_per_Cu, preserved), where list_per_Cu is a list of 6-distance lists.
    Octahedral preserved if for every Cu, 6th-nearest O distance ≤ 3.0 Å and
    a clear 4+2 split exists: d5 - d4 ≥ 0.05 Å.
    """
    symbols = atoms.get_chemical_symbols()
    cu_indices = [i for i, s in enumerate(symbols) if s == "Cu"]
    o_indices = [i for i, s in enumerate(symbols) if s == "O"]

    all_dists: List[List[float]] = []
    preserved = True
    for i in cu_indices:
        dists = [atoms.get_distance(i, j, mic=True) for j in o_indices]
        dists.sort()
        six = dists[:6] if len(dists) >= 6 else dists
        if len(six) < 6:
            preserved = False
        else:
            # Check thresholds
            if six[5] > 3.0:
                preserved = False
            if (six[4] - six[3]) < 0.05:
                preserved = False
        all_dists.append(six)
    return all_dists, preserved


def aggregate_bonds(cu_o_lists: List[List[float]]) -> Dict[str, float]:
    """Aggregate Cu–O distances into equatorial (first 4) and apical (last 2).

    Returns dict with mean and std for equatorial and apical groups.
    """
    eq = []
    ap = []
    for six in cu_o_lists:
        if len(six) >= 6:
            eq.extend(six[:4])
            ap.extend(six[4:6])
    res = {}
    if eq:
        res["eq_mean_A"] = float(np.mean(eq))
        res["eq_std_A"] = float(np.std(eq))
    else:
        res["eq_mean_A"] = float("nan")
        res["eq_std_A"] = float("nan")
    if ap:
        res["ap_mean_A"] = float(np.mean(ap))
        res["ap_std_A"] = float(np.std(ap))
    else:
        res["ap_mean_A"] = float("nan")
        res["ap_std_A"] = float("nan")
    return res


@dataclass
class Result:
    c_scale: float
    a_A: float
    b_A: float
    c_A: float
    natoms: int
    energy_eV: float
    energy_per_atom_eV: float
    max_force_eV_per_A: float
    octahedral_preserved: bool
    eq_mean_A: float
    eq_std_A: float
    ap_mean_A: float
    ap_std_A: float
    cu_o_bonds_A: List[List[float]]


def main():
    # Lazy imports for heavy libs; load CHGNet once
    from pymatgen.io.ase import AseAtomsAdaptor  # type: ignore
    from ase.io import write as ase_write
    from chgnet.model import CHGNet, CHGNetCalculator  # type: ignore

    # 1) Get starting structure
    pmg_struct = get_initial_structure()

    # 2) Convert to ASE Atoms
    adaptor = AseAtomsAdaptor()
    atoms0 = adaptor.get_atoms(pmg_struct)

    # 3) Load CHGNet once and reuse calculator
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    expansions = [1.00, 1.02, 1.05, 1.10]
    results: List[Result] = []

    # Record baseline lattice parameters
    base_cell = atoms0.cell.lengths()

    # Loop over c-axis scale factors
    for s in expansions:
        atoms = atoms_copy_with_c_scale(atoms0, s)
        atoms.set_calculator(calc)

        # Guard against NaNs before optimization
        ensure_finite_energy_forces(atoms, tag=f"before_opt_s={s}")

        # Relax positions only (two-stage BFGS)
        relax_positions_only(atoms, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300, log_prefix=f"relax_c{s:.3f}")

        # Guard after optimization
        ensure_finite_energy_forces(atoms, tag=f"after_opt_s={s}")

        # Compute properties
        e = float(atoms.get_potential_energy())
        n = len(atoms)
        epa = e / n
        fmax = compute_max_force(atoms)
        cuo_lists, preserved = cu_o_bonds(atoms)
        agg = aggregate_bonds(cuo_lists)

        # Lattice parameters
        aA, bA, cA = atoms.cell.lengths()

        # Save final structure
        poscar_name = f"POSCAR_La2CuO4_cscale_{s:.3f}.vasp"
        ase_write(poscar_name, atoms, format="vasp", direct=True, vasp5=True)

        res = Result(
            c_scale=s,
            a_A=float(aA),
            b_A=float(bA),
            c_A=float(cA),
            natoms=n,
            energy_eV=e,
            energy_per_atom_eV=float(epa),
            max_force_eV_per_A=float(fmax),
            octahedral_preserved=bool(preserved),
            eq_mean_A=float(agg["eq_mean_A"]),
            eq_std_A=float(agg["eq_std_A"]),
            ap_mean_A=float(agg["ap_mean_A"]),
            ap_std_A=float(agg["ap_std_A"]),
            cu_o_bonds_A=[[float(x) for x in six] for six in cuo_lists],
        )
        results.append(res)

    # 4) Write results.json and summary.csv
    out_json = "results_la2cuo4_cscan.json"
    with open(out_json, "w") as f:
        json.dump([asdict(r) for r in results], f, indent=2)

    # CSV summary
    import csv

    out_csv = "summary_la2cuo4_cscan.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "c_scale",
            "a_A",
            "b_A",
            "c_A",
            "natoms",
            "energy_eV",
            "energy_per_atom_eV",
            "max_force_eV_per_A",
            "octahedral_preserved",
            "eq_mean_A",
            "eq_std_A",
            "ap_mean_A",
            "ap_std_A",
        ])
        for r in results:
            w.writerow([
                r.c_scale,
                r.a_A,
                r.b_A,
                r.c_A,
                r.natoms,
                r.energy_eV,
                r.energy_per_atom_eV,
                r.max_force_eV_per_A,
                r.octahedral_preserved,
                r.eq_mean_A,
                r.eq_std_A,
                r.ap_mean_A,
                r.ap_std_A,
            ])

    # 5) Print concise human summary
    print("La2CuO4 c-axis scan (positions-only CHGNet relax):")
    print("c_scale  a(Å)  c(Å)  E/atom(eV)  Fmax(eV/Å)  eq_mean(Å)  ap_mean(Å)  oct_preserved")
    for r in results:
        print(
            f"{r.c_scale:>6.3f}  {r.a_A:>4.2f}  {r.c_A:>5.2f}  {r.energy_per_atom_eV:>10.5f}  "
            f"{r.max_force_eV_per_A:>10.4f}  {r.eq_mean_A:>9.3f}  {r.ap_mean_A:>9.3f}  {str(r.octahedral_preserved):>5}"
        )


if __name__ == "__main__":
    main()

