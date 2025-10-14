#!/usr/bin/env python3
import os
import json
import csv
import time
import random
import numpy as np

from ase.build import bulk, make_supercell
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator


# Units
EV_PER_A3_TO_GPA = 160.21766208


def set_seeds(seed: int = 42):
    random.seed(seed)
    np.random.seed(seed)


def check_finite(atoms):
    import numpy as _np
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not _np.isfinite(e) or not _np.isfinite(f).all():
        raise RuntimeError("Non-finite energy/forces; check initial structure or calculator state.")


def relax_positions(atoms, calc, fmax=0.1, steps=150, logfile=None):
    atoms.set_calculator(calc)
    check_finite(atoms)
    opt = BFGS(atoms, logfile=logfile)
    opt.run(fmax=fmax, steps=steps)
    # Final sanity check
    check_finite(atoms)
    return atoms


def relax_cell_hydrostatic(atoms, calc, fmax=0.05, steps=300, logfile=None):
    atoms.set_calculator(calc)
    check_finite(atoms)
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt = BFGS(ucf, logfile=logfile)
    opt.run(fmax=fmax, steps=steps)
    check_finite(atoms)
    return atoms


def voigt_to_strain_matrix(e6):
    """Map 6-component Voigt strain vector to 3x3 symmetric small-strain tensor.
    Voigt convention: e = [e11, e22, e33, e23(=gamma23), e13, e12], with shear components being engineering strains.
    So epsilon_23 = e4/2, etc.
    """
    e11, e22, e33, e23, e13, e12 = e6
    eps = np.array(
        [
            [e11, e12 / 2.0, e13 / 2.0],
            [e12 / 2.0, e22, e23 / 2.0],
            [e13 / 2.0, e23 / 2.0, e33],
        ]
    )
    return eps


def apply_strain(atoms, eps):
    new = atoms.copy()
    F = np.eye(3) + eps  # small-strain approximation
    new.set_cell(F @ new.cell, scale_atoms=True)
    return new


def stress_voigt_GPa(atoms):
    # ASE returns stress as [sxx, syy, szz, syz, sxz, sxy] in eV/Å^3 (sign convention: positive for compressive)
    s = atoms.get_stress(voigt=True)
    return np.array(s) * EV_PER_A3_TO_GPA


def stress_for_F(atoms_base, calc, F, fmax=0.10, steps=150):
    a = atoms_base.copy()
    a.set_cell(F @ a.cell, scale_atoms=True)
    a = relax_positions(a, calc, fmax=fmax, steps=steps, logfile=None)
    return stress_voigt_GPa(a)


def energy_density_GPa(atoms, V0=None):
    e = atoms.get_potential_energy()  # eV
    if V0 is None:
        V0 = atoms.get_volume()
    U = (e / V0) * EV_PER_A3_TO_GPA    # GPa, per reference (unstrained) volume
    return U


def compute_cubic_moduli_by_energy(atoms_relaxed, calc, deltas=(0.004, 0.006), per_strain_relax_fmax=0.10, per_strain_relax_steps=150):
    """Compute effective cubic elastic constants via energy second differences.
    Returns list of dicts per delta: {delta, K, C11, C12, C44, G_VRH, G_RH, G, E, nu}.
    """
    base = atoms_relaxed.copy()
    base.set_calculator(calc)
    V0 = base.get_volume()
    U0 = energy_density_GPa(base, V0=V0)

    out = []
    I = np.eye(3)
    for d in deltas:
        # Hydrostatic: F = (1+δ) I
        Fp = (1.0 + d) * I
        Fm = (1.0 - d) * I
        ap = base.copy(); ap.set_cell(Fp @ ap.cell, scale_atoms=True); ap.set_calculator(calc)
        am = base.copy(); am.set_cell(Fm @ am.cell, scale_atoms=True); am.set_calculator(calc)
        Up = energy_density_GPa(ap, V0=V0)
        Um = energy_density_GPa(am, V0=V0)
        # K from U(δ) = (9/2) K δ^2; using central difference coefficient b = (U+ + U- - 2U0) / (2 δ^2)
        b_h = (Up + Um - 2.0 * U0) / (2.0 * d * d)
        K = (2.0 / 9.0) * b_h

        # Orthorhombic: e1=+δ, e2=-δ, e3=0  => F = diag(1+δ, 1-δ, 1)
        Fo_p = np.diag([1.0 + d, 1.0 - d, 1.0])
        Fo_m = np.diag([1.0 - d, 1.0 + d, 1.0])
        aop = base.copy(); aop.set_cell(Fo_p @ aop.cell, scale_atoms=True); aop.set_calculator(calc)
        aom = base.copy(); aom.set_cell(Fo_m @ aom.cell, scale_atoms=True); aom.set_calculator(calc)
        Uop = energy_density_GPa(aop, V0=V0)
        Uom = energy_density_GPa(aom, V0=V0)
        b_o = (Uop + Uom - 2.0 * U0) / (2.0 * d * d)
        A = b_o  # A = C11 - C12

        # Monoclinic shear (pure symmetric shear): eps12 = d/2
        eps = np.zeros((3, 3)); eps[0, 1] = d / 2.0; eps[1, 0] = d / 2.0
        Fs_p = I + eps
        Fs_m = I - eps
        asp = base.copy(); asp.set_cell(Fs_p @ asp.cell, scale_atoms=True); asp.set_calculator(calc)
        asm = base.copy(); asm.set_cell(Fs_m @ asm.cell, scale_atoms=True); asm.set_calculator(calc)
        Usp = energy_density_GPa(asp, V0=V0)
        Usm = energy_density_GPa(asm, V0=V0)
        b_s = (Usp + Usm - 2.0 * U0) / (2.0 * d * d)
        C44 = 2.0 * b_s

        # Recover C11 and C12 from K and A
        C11 = (2.0 * A + 3.0 * K) / 3.0
        C12 = (3.0 * K - A) / 3.0

        # VRH shear for cubic
        G_V = (C11 - C12 + 3.0 * C44) / 5.0
        G_R = 5.0 * (C11 - C12) * C44 / (4.0 * C44 + 3.0 * (C11 - C12))
        G = 0.5 * (G_V + G_R)
        E = 9.0 * K * G / (3.0 * K + G)
        nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))

        out.append({
            "delta": d,
            "K": float(K),
            "C11": float(C11),
            "C12": float(C12),
            "C44": float(C44),
            "G_V": float(G_V),
            "G_R": float(G_R),
            "G": float(G),
            "E": float(E),
            "nu": float(nu),
        })
    return out


def vrh_moduli_from_C(C):
    """Compute Voigt-Reuss-Hill bulk, shear, Young's modulus and Poisson's ratio from a 6x6 stiffness matrix C (GPa)."""
    # Ensure array and symmetry
    C = 0.5 * (C + C.T)
    # Compliance
    try:
        S = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        # Add small regularization if nearly singular
        eigvals, eigvecs = np.linalg.eigh(C)
        eigvals_reg = np.clip(eigvals, 1e-6, None)
        C_reg = (eigvecs @ np.diag(eigvals_reg) @ eigvecs.T)
        S = np.linalg.inv(C_reg)

    # Voigt averages
    K_V = (C[0, 0] + C[1, 1] + C[2, 2] + 2.0 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9.0
    G_V = (
        (C[0, 0] + C[1, 1] + C[2, 2])
        - (C[0, 1] + C[0, 2] + C[1, 2])
        + 3.0 * (C[3, 3] + C[4, 4] + C[5, 5])
    ) / 15.0

    # Reuss averages
    denom_K = S[0, 0] + S[1, 1] + S[2, 2] + 2.0 * (S[0, 1] + S[0, 2] + S[1, 2])
    K_R = 1.0 / denom_K
    denom_G = 4.0 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4.0 * (S[0, 1] + S[0, 2] + S[1, 2]) + 3.0 * (S[3, 3] + S[4, 4] + S[5, 5])
    G_R = 15.0 / denom_G

    K = 0.5 * (K_V + K_R)
    G = 0.5 * (G_V + G_R)

    E = 9.0 * K * G / (3.0 * K + G)
    nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G))
    return {
        "K_V": float(K_V),
        "G_V": float(G_V),
        "K_R": float(K_R),
        "G_R": float(G_R),
        "K": float(K),
        "G": float(G),
        "E": float(E),
        "nu": float(nu),
    }


def build_pure_al(a0=4.05, super=(3, 3, 3)):
    # Conventional cubic fcc cell (4 atoms), then supercell
    al_conv = bulk("Al", crystalstructure="fcc", a=a0, cubic=True)
    P = np.diag(super)
    al_sc = make_supercell(al_conv, P)
    return al_sc


def build_al_mg_solution(a0=4.05, super=(3, 3, 3), n_mg=2, seed=42):
    atoms = build_pure_al(a0=a0, super=super)
    set_seeds(seed)
    N = len(atoms)
    assert n_mg < N
    idxs = np.random.choice(np.arange(N), size=n_mg, replace=False)
    symbols = atoms.get_chemical_symbols()
    for i in idxs:
        symbols[i] = "Mg"
    atoms.set_chemical_symbols(symbols)
    return atoms


def main():
    set_seeds(42)

    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    out_json = "results.json"
    out_csv = "summary.csv"

    systems = []

    # System 1: Pure Al
    al = build_pure_al(a0=4.05, super=(1, 1, 1))  # 4 atoms (cubic conventional cell)
    al.set_calculator(calc)

    # Stage 1: positions-only relax
    al = relax_positions(al, calc, fmax=0.2, steps=150, logfile="al_relax1.log")
    # Stage 2: hydrostatic cell relax
    al = relax_cell_hydrostatic(al, calc, fmax=0.05, steps=300, logfile="al_relax2.log")
    write("Al_relaxed_POSCAR", al, format="vasp")

    # Elastic moduli via FD (two deltas for an uncertainty estimate)
    al_mod = compute_cubic_moduli_by_energy(al, calc, deltas=(0.01, 0.015), per_strain_relax_fmax=0.10, per_strain_relax_steps=150)
    al_Es = [d["E"] for d in al_mod]
    al_result = {
        "name": "Al_fcc",
        "formula": al.get_chemical_formula(),
        "n_atoms": len(al),
        "moduli": al_mod,
        "aggregate": {
            "E_mean_GPa": float(np.mean(al_Es)),
            "E_std_GPa": float(np.std(al_Es, ddof=1) if len(al_Es) > 1 else 0.0),
            "K_mean_GPa": float(np.mean([d["K"] for d in al_mod])),
            "G_mean_GPa": float(np.mean([d["G"] for d in al_mod])),
            "nu_mean": float(np.mean([d["nu"] for d in al_mod])),
        },
    }
    systems.append(al_result)

    # System 2: Al-Mg ~1.85 at% (2 Mg in 108)
    almg = build_al_mg_solution(a0=4.05, super=(3, 3, 3), n_mg=2, seed=42)
    almg.set_calculator(calc)
    # Stage 1: positions-only relax
    almg = relax_positions(almg, calc, fmax=0.2, steps=200, logfile="almg_relax1.log")
    # Stage 2: hydrostatic cell relax
    almg = relax_cell_hydrostatic(almg, calc, fmax=0.05, steps=350, logfile="almg_relax2.log")
    write("AlMg_relaxed_POSCAR", almg, format="vasp")

    # Elastic moduli (two deltas)
    almg_mod = compute_cubic_moduli_by_energy(almg, calc, deltas=(0.01, 0.015), per_strain_relax_fmax=0.10, per_strain_relax_steps=160)
    almg_Es = [d["E"] for d in almg_mod]
    almg_result = {
        "name": "Al-1.85%Mg_fcc",
        "formula": almg.get_chemical_formula(),
        "n_atoms": len(almg),
        "moduli": almg_mod,
        "aggregate": {
            "E_mean_GPa": float(np.mean(almg_Es)),
            "E_std_GPa": float(np.std(almg_Es, ddof=1) if len(almg_Es) > 1 else 0.0),
            "K_mean_GPa": float(np.mean([d["K"] for d in almg_mod])),
            "G_mean_GPa": float(np.mean([d["G"] for d in almg_mod])),
            "nu_mean": float(np.mean([d["nu"] for d in almg_mod])),
        },
    }
    systems.append(almg_result)

    # Write summary CSV (per-delta rows)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["system", "delta", "K(GPa)", "C11(GPa)", "C12(GPa)", "C44(GPa)", "G(GPa)", "E(GPa)", "nu"])
        for sys in systems:
            for row in sys["moduli"]:
                w.writerow([sys["name"], row["delta"], row["K"], row.get("C11"), row.get("C12"), row.get("C44"), row["G"], row["E"], row["nu"]])

    # Write JSON with details
    payload = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "systems": systems,
    }
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=2)

    # Human summary
    E_al = al_result["aggregate"]["E_mean_GPa"]
    dE_al = al_result["aggregate"]["E_std_GPa"]
    E_almg = almg_result["aggregate"]["E_mean_GPa"]
    dE_almg = almg_result["aggregate"]["E_std_GPa"]
    change_pct = 100.0 * (E_almg - E_al) / E_al

    print("\nResults (VRH averages):")
    print(f"- Al (fcc, {al_result['n_atoms']} atoms): E = {E_al:.2f} +/- {dE_al:.2f} GPa")
    print(f"- Al-1.85%Mg (fcc, {almg_result['n_atoms']} atoms): E = {E_almg:.2f} +/- {dE_almg:.2f} GPa")
    print(f"- Relative change: {change_pct:.2f}%")

    if change_pct >= 5.0:
        print("Conclusion: Dissolved ~2 at% Mg significantly increases (>= 5%) the Young's modulus.")
    else:
        print("Conclusion: Dissolved ~2 at% Mg does NOT increase Young's modulus by >= 5%.")


if __name__ == "__main__":
    main()
