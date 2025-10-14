#!/usr/bin/env python3
"""
CHGNet single-defect relaxation in fcc Al and lattice parameter comparison.

Workflow (single script, single model load):
  1) Build fcc Al conventional cell, create 2x2x2 supercell (32 atoms)
  2) For each case (pure, Cu, Mg, Zn):
       - Substitute 1 Al with dopant (pure keeps Al)
       - Stage 1 relax: atomic positions only (BFGS)
       - Stage 2 relax: isotropic cell (UnitCellFilter with hydrostatic_strain=True)
       - Record final lattice parameter a = mean(cell lengths)/rep
  3) Write final POSCARs, results.json, summary.csv
  4) Print concise human summary and Δa vs pure baseline

Notes:
  - Model loaded once and reused across all calculations.
  - Conservative settings for robustness and speed; isotropic cell relax gives a well-defined a.
  - Seeded for determinism.
"""

import os
import json
import csv
import math
import time
import numpy as np

# Determinism
np.random.seed(42)
try:
    import torch
    torch.manual_seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(42)
except Exception:
    torch = None  # torch not strictly required for running CHGNetCalculator

from ase.build import bulk
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator


def cartesian_center(cell):
    # Center of the cell in Cartesian coords: M*[0.5,0.5,0.5]
    mat = cell.array  # 3x3
    return mat.dot(np.array([0.5, 0.5, 0.5]))


def nearest_atom_index_to_center(atoms):
    c = cartesian_center(atoms.get_cell())
    pos = atoms.get_positions()
    d2 = np.sum((pos - c) ** 2, axis=1)
    return int(np.argmin(d2))


def max_force(atoms):
    f = atoms.get_forces()
    return float(np.max(np.linalg.norm(f, axis=1)))


def ensure_finite(atoms, where=""):
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not np.isfinite(e) or not np.isfinite(f).all():
        raise RuntimeError(f"Non-finite energy/forces {where}")
    return float(e), float(np.max(np.linalg.norm(f, axis=1)))


def lattice_parameter_from_cell(atoms, reps):
    # Effective cubic parameter from the (possibly slightly distorted) cell.
    # Use mean of the three lengths, divided by repetition along one axis.
    L = atoms.cell.lengths()
    a_eff = float(np.mean(L) / reps[0])
    return a_eff


def relax_structure(atoms, calc, label, reps=(2, 2, 2),
                    fmax_coarse=0.2, fmax_fine=0.05,
                    steps1=150, steps2=200, workdir="."):
    """Two-stage relaxation: positions-only -> hydrostatic cell.

    Returns dict with energies, fmax, steps, a_final, stress, etc.
    """
    atoms.set_calculator(calc)

    # Stage 0 checks
    e0, fmax0 = ensure_finite(atoms, where="(pre-opt)")

    # Stage 1: positions-only
    log1 = os.path.join(workdir, f"{label}_relax1.log")
    opt1 = BFGS(atoms, logfile=log1)
    t1_start = time.time()
    opt1.run(fmax=fmax_coarse, steps=steps1)
    t1 = time.time() - t1_start
    e1, fmax1 = ensure_finite(atoms, where="(post opt1)")
    nsteps1 = int(getattr(opt1, 'nsteps', -1))

    # Stage 2: hydrostatic cell relax
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    log2 = os.path.join(workdir, f"{label}_relax2.log")
    opt2 = BFGS(ucf, logfile=log2)
    t2_start = time.time()
    opt2.run(fmax=fmax_fine, steps=steps2)
    t2 = time.time() - t2_start
    e2, fmax2 = ensure_finite(atoms, where="(post opt2)")
    nsteps2 = int(getattr(opt2, 'nsteps', -1))

    # Final descriptors
    a_final = lattice_parameter_from_cell(atoms, reps)
    try:
        stress = atoms.get_stress(voigt=True)
        stress = [float(x) for x in stress]
    except Exception:
        stress = None

    return {
        "label": label,
        "natoms": int(len(atoms)),
        "a_final": a_final,
        "energies": {
            "initial": e0,
            "after_stage1": e1,
            "after_stage2": e2,
            "delta1": float(e1 - e0),
            "delta2": float(e2 - e1),
            "total_drop": float(e2 - e0),
        },
        "fmax": {
            "initial": fmax0,
            "after_stage1": fmax1,
            "after_stage2": fmax2,
        },
        "steps": {"stage1": nsteps1, "stage2": nsteps2},
        "timing_s": {"stage1": t1, "stage2": t2},
        "stress_voigt_final": stress,
    }


def main():
    workdir = os.getcwd()
    reps = (2, 2, 2)  # 32-atom supercell, ~3.125 at% single substitution
    a_init = 4.05  # Al fcc starting guess (Å)

    # Build starting supercell once
    base = bulk('Al', 'fcc', a=a_init, cubic=True)
    pure = base.repeat(reps)

    # Identify a reproducible substitution site (near geometric center)
    sub_index = nearest_atom_index_to_center(pure)

    # Prepare structures
    systems = {
        "pure": pure.copy(),
        "Cu": pure.copy(),
        "Mg": pure.copy(),
        "Zn": pure.copy(),
    }
    for dop in ["Cu", "Mg", "Zn"]:
        sy = systems[dop]
        sy[sub_index].symbol = dop

    # Load CHGNet once, reuse calculator
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    results = []
    final_structures = {}

    # Relax pure first to establish baseline a
    order = ["pure", "Cu", "Mg", "Zn"]

    for label in order:
        atoms = systems[label]

        # Clip extremely close pairs if any (safety; scale up slightly)
        # For fcc + single substitution, this should not trigger.
        # Min distance threshold (Å)
        mind = 0.7
        pos = atoms.get_positions()
        # Simple O(N^2) check for small systems
        bad = False
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                rij = np.linalg.norm(pos[i] - pos[j])
                if rij < mind:
                    bad = True
                    break
            if bad:
                break
        if bad:
            atoms.set_cell(atoms.get_cell() * 1.05, scale_atoms=True)

        res = relax_structure(atoms, calc, label, reps=reps, workdir=workdir)
        results.append(res)
        final_structures[label] = atoms.copy()

        # Write POSCAR for each relaxed system
        poscar_name = f"POSCAR_Al_{label}.vasp" if label != "pure" else "POSCAR_Al_pure.vasp"
        write(os.path.join(workdir, poscar_name), atoms, format='vasp', direct=True, vasp5=True)

    # Compute Δa relative to pure
    a_pure = next(r["a_final"] for r in results if r["label"] == "pure")
    for r in results:
        da = r["a_final"] - a_pure
        pct = (da / a_pure) * 100.0
        r["delta_a_A"] = float(da)
        r["delta_percent"] = float(pct)

    # Write results.json
    with open(os.path.join(workdir, "results.json"), "w") as f:
        json.dump({"reps": reps, "a_init": a_init, "sub_index": sub_index, "results": results}, f, indent=2)

    # Write summary.csv
    with open(os.path.join(workdir, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "natoms", "a_final_A", "delta_a_A", "delta_percent", "steps1", "steps2"]) 
        for r in results:
            w.writerow([
                r["label"], r["natoms"], f"{r['a_final']:.6f}", f"{r['delta_a_A']:.6f}", f"{r['delta_percent']:.4f}",
                r["steps"]["stage1"], r["steps"]["stage2"],
            ])

    # Human-readable summary
    print("\n=== CHGNet Al fcc single-substitution relaxations (2x2x2 supercell) ===")
    print(f"Baseline (pure Al) lattice parameter a = {a_pure:.6f} Å")
    print("\nFinal lattice parameters and changes:")
    print(f"{'Label':<6} {'a (Å)':>12} {'Δa (Å)':>12} {'Δa (%)':>10} {'Steps1':>8} {'Steps2':>8}")
    for r in results:
        print(f"{r['label']:<6} {r['a_final']:>12.6f} {r['delta_a_A']:>12.6f} {r['delta_percent']:>10.4f} {r['steps']['stage1']:>8d} {r['steps']['stage2']:>8d}")

    # Rank dopants by magnitude of lattice expansion (absolute Δa), exclude pure
    dopants_only = [r for r in results if r["label"] != "pure"]
    ranked = sorted(dopants_only, key=lambda x: abs(x["delta_a_A"]), reverse=True)
    print("\nRanking by |Δa| (largest expansion magnitude first):")
    for i, r in enumerate(ranked, 1):
        print(f"{i}. {r['label']}: Δa = {r['delta_a_A']:.6f} Å ({r['delta_percent']:.4f}%)")

    # Also print sign-based ranking (expansion > contraction), largest positive first
    ranked_pos = sorted(dopants_only, key=lambda x: x["delta_a_A"], reverse=True)
    print("\nRanking by Δa (expansion → contraction):")
    for i, r in enumerate(ranked_pos, 1):
        print(f"{i}. {r['label']}: Δa = {r['delta_a_A']:.6f} Å ({r['delta_percent']:.4f}%)")


if __name__ == "__main__":
    main()

