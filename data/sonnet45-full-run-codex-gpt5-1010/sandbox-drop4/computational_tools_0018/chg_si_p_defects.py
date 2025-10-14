import os
import json
import csv
import math
import numpy as np

from ase.build import bulk
from ase import Atoms
from ase.atom import Atom
from ase.optimize import BFGS

from chgnet.model import CHGNet, CHGNetCalculator


def max_force_magnitude(atoms: Atoms) -> float:
    f = atoms.get_forces()
    if not np.isfinite(f).all():
        return np.inf
    return float(np.max(np.linalg.norm(f, axis=1)))


def ensure_finite(atoms: Atoms, context: str = ""):
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not (np.isfinite(e) and np.isfinite(f).all()):
        raise RuntimeError(f"Non-finite energy/forces {context}.")


def relax_positions(atoms: Atoms, calc: CHGNetCalculator, label: str,
                    fmax_coarse: float = 0.2, fmax_fine: float = 0.05,
                    steps_coarse: int = 75, steps_fine: int = 150):
    atoms.set_calculator(calc)
    atoms.set_pbc((True, True, True))

    fmax_log = []

    def log_fmax():
        try:
            fmax_log.append(max_force_magnitude(atoms))
        except Exception:
            # Fail-safe: don't break optimization if logging fails
            pass

    # Pre-check
    ensure_finite(atoms, context=f"before relax ({label})")

    # Stage 1: positions-only, coarse
    opt1 = BFGS(atoms, logfile=f"{label}_relax1.log")
    opt1.attach(log_fmax, interval=1)
    opt1.run(fmax=fmax_coarse, steps=steps_coarse)

    ensure_finite(atoms, context=f"after stage 1 ({label})")

    # Stage 2: positions-only, fine
    opt2 = BFGS(atoms, logfile=f"{label}_relax2.log")
    opt2.attach(log_fmax, interval=1)
    opt2.run(fmax=fmax_fine, steps=steps_fine)

    ensure_finite(atoms, context=f"after stage 2 ({label})")

    e = float(atoms.get_potential_energy())
    fmax_final = max_force_magnitude(atoms)

    return {
        "final_energy": e,
        "fmax_trajectory": fmax_log,
        "final_fmax": float(fmax_final),
        "label": label,
    }


def build_si_supercell(a0: float = 5.431, reps=(2, 2, 2)) -> Atoms:
    # Conventional cubic diamond cell (8 atoms) then 2x2x2 -> 64 atoms
    si = bulk("Si", crystalstructure="diamond", a=a0, cubic=True)
    si_sc = si.repeat(reps)
    si_sc.set_pbc((True, True, True))
    return si_sc


def make_substitutional_P(si_sc: Atoms) -> Atoms:
    sub = si_sc.copy()
    # Replace one Si with P (index 0 for simplicity)
    sub[0].symbol = "P"
    return sub


def frac_to_cart(cell: np.ndarray, frac: np.ndarray) -> np.ndarray:
    # frac is shape (3,), cell is 3x3
    return np.dot(frac, cell)


def make_interstitial_P(si_sc: Atoms) -> Atoms:
    inter = si_sc.copy()
    cell = np.array(inter.cell)  # 3x3
    # Tetrahedral interstitial in diamond: 8a Wyckoff ~ (1/8,1/8,1/8) of the conventional cell.
    # For a 2x2x2 supercell, a corresponding site near the supercell center is (0.5+1/16) along each axis.
    frac = np.array([0.5625, 0.5625, 0.5625])  # 0.5 + 1/16
    pos = frac_to_cart(cell, frac)
    inter.append(Atom("P", position=pos))
    return inter


def write_structures(bulk_sc: Atoms, sub: Atoms, inter: Atoms):
    from ase.io import write
    write("POSCAR_bulk.vasp", bulk_sc, format="vasp")
    write("POSCAR_sub.vasp", sub, format="vasp")
    write("POSCAR_int.vasp", inter, format="vasp")


def main():
    # Reproducibility: set seeds
    try:
        import torch
        torch.manual_seed(42)
    except Exception:
        pass
    np.random.seed(42)

    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Build structures
    bulk_sc = build_si_supercell(a0=5.431, reps=(2, 2, 2))  # 64 atoms
    sub = make_substitutional_P(bulk_sc)
    inter = make_interstitial_P(bulk_sc)

    # Relax structures (positions only, consistent settings)
    res_bulk = relax_positions(bulk_sc, calc, label="bulk")
    res_sub = relax_positions(sub, calc, label="sub")
    res_int = relax_positions(inter, calc, label="int")

    # Compute energies and ΔE
    e_bulk = res_bulk["final_energy"]
    e_sub = res_sub["final_energy"]
    e_int = res_int["final_energy"]

    n_bulk = len(bulk_sc)  # 64
    mu_si = e_bulk / n_bulk  # eV/atom

    # ΔE = E_f(P_i) - E_f(P_Si) = E_tot(Si_N+P) - E_tot(Si_{N-1}P) - μ_Si
    delta_E = e_int - e_sub - mu_si

    # Write structures
    write_structures(bulk_sc, sub, inter)

    # Prepare results
    results = {
        "method": "CHGNet",
        "supercell": "Si diamond conventional 2x2x2 (64 atoms)",
        "a0_Si_A": 5.431,
        "energies_eV": {
            "E_bulk_64": e_bulk,
            "E_substitutional": e_sub,
            "E_interstitial": e_int,
            "mu_Si_per_atom": mu_si,
        },
        "delta_E_eV": delta_E,
        "notes": "ΔE = E_f(P_i) - E_f(P_Si)",
        "relaxation": {
            "bulk": {"final_fmax": res_bulk["final_fmax"], "steps": len(res_bulk["fmax_trajectory"])},
            "sub": {"final_fmax": res_sub["final_fmax"], "steps": len(res_sub["fmax_trajectory"])},
            "int": {"final_fmax": res_int["final_fmax"], "steps": len(res_int["fmax_trajectory"])},
        },
    }

    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)

    # CSV summary
    with open("summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value_eV"])
        w.writerow(["E_bulk_64", e_bulk])
        w.writerow(["mu_Si_per_atom", mu_si])
        w.writerow(["E_substitutional", e_sub])
        w.writerow(["E_interstitial", e_int])
        w.writerow(["delta_E = Ef(Pi)-Ef(Psi)", delta_E])

    # Human-readable concise summary
    print("Computed with CHGNet on Si 2x2x2 supercell (64 atoms).")
    print(f"E_bulk(64) = {e_bulk:.6f} eV; mu_Si = {mu_si:.6f} eV/atom")
    print(f"E_sub = {e_sub:.6f} eV; E_int = {e_int:.6f} eV")
    print(f"ΔE = E_f(P_i) - E_f(P_Si) = {delta_E:.6f} eV")


if __name__ == "__main__":
    main()

