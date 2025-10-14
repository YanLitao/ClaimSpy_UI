import os
import json
import numpy as np
from dataclasses import asdict, dataclass
from collections import Counter

# Ensure we do not accidentally source user profiles that slow runs
os.environ.setdefault("BASH_ENV", "/dev/null")

# Optional: load MP_API_KEY if present (not used here but allowed by instructions)
try:
    from dotenv import load_dotenv
    load_dotenv(".env")
except Exception:
    pass

from ase.build import bulk
from ase.io import write
from ase import Atoms
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.units import kJ

from chgnet.model import CHGNet, CHGNetCalculator


GPA_PER_eV_per_A3 = 160.21766208


def seed_all(seed: int = 42):
    np.random.seed(seed)


def make_fcc_al(a0: float = 4.05) -> Atoms:
    at = bulk("Al", "fcc", a=a0, cubic=True)
    at.pbc = True
    return at


def supercell(at: Atoms, reps=(2, 2, 2)) -> Atoms:
    sc = at.repeat(reps)
    sc.pbc = True
    return sc


def substitute_cu(at: Atoms, n_cu: int, min_sep: float = 7.5) -> Atoms:
    """Substitute n_cu Al atoms with Cu, enforcing a minimum pair separation (Å)."""
    positions = at.get_positions()
    N = len(at)
    chosen = []
    # Greedy farthest-point sampling
    for _ in range(n_cu):
        candidates = np.setdiff1d(np.arange(N), np.array(chosen, dtype=int))
        np.random.shuffle(candidates)
        placed = False
        for idx in candidates:
            pos = positions[idx]
            if all(np.linalg.norm(pos - positions[j]) >= min_sep for j in chosen):
                chosen.append(int(idx))
                placed = True
                break
        if not placed:
            # fallback: relax criterion slightly
            for idx in candidates:
                pos = positions[idx]
                if all(np.linalg.norm(pos - positions[j]) >= (0.5 * min_sep) for j in chosen):
                    chosen.append(int(idx))
                    placed = True
                    break
        if not placed:
            # last resort: just pick the farthest from the existing set
            if not chosen:
                chosen.append(int(candidates[0]))
            else:
                dists = []
                for idx in candidates:
                    pos = positions[idx]
                    d = min(np.linalg.norm(pos - positions[j]) for j in chosen)
                    dists.append((d, int(idx)))
                dists.sort(reverse=True)
                chosen.append(dists[0][1])
    symbols = at.get_chemical_symbols()
    for idx in chosen:
        symbols[idx] = "Cu"
    at.set_chemical_symbols(symbols)
    return at


def ensure_min_distances(at: Atoms, rmin: float = 0.7):
    """If any pair is closer than rmin Å, isotropically expand the cell and positions."""
    from ase.neighborlist import neighbor_list
    scale = 1.0
    for _ in range(5):
        i, j, d = neighbor_list("ijd", at, cutoff=rmin)
        if len(d) == 0:
            return
        # Expand
        scale *= 1.05
        cell = at.get_cell()
        at.set_cell(cell * 1.05, scale_atoms=True)
    # Final check; if still problematic, raise
    i, j, d = neighbor_list("ijd", at, cutoff=rmin)
    if len(d) != 0:
        raise RuntimeError("Geometry contains unphysically short distances even after expansion.")


def relax_structure(at: Atoms, calc: CHGNetCalculator, tag: str,
                    fmax_coarse: float = 0.2, fmax_fine: float = 0.05,
                    steps_coarse: int = 150, steps_fine: int = 300,
                    variable_cell: bool = True) -> Atoms:
    at = at.copy()
    at.calc = calc
    # Stage 1: positions only
    opt1 = BFGS(at, logfile=f"relax1_{tag}.log")
    opt1.run(fmax=fmax_coarse, steps=steps_coarse)
    # Stage 2: variable cell (hydrostatic) if requested
    if variable_cell:
        # Allow full cell degrees of freedom to reach a true zero-stress state
        ucf = UnitCellFilter(at)
        opt2 = BFGS(ucf, logfile=f"relax2_{tag}.log")
        opt2.run(fmax=fmax_fine, steps=steps_fine)
    else:
        opt2 = BFGS(at, logfile=f"relax2_{tag}.log")
        opt2.run(fmax=fmax_fine, steps=steps_fine)
    return at


def strain_matrix_from_voigt(eps: np.ndarray) -> np.ndarray:
    """Convert 6-strain Voigt to 3x3 symmetric small-strain tensor.
    Voigt order: [e_xx, e_yy, e_zz, e_yz, e_xz, e_xy]
    with engineering shear components: e_yz = 2*epsilon_yz, etc.
    """
    exx, eyy, ezz, eyz, exz, exy = eps
    E = np.array([
        [exx, 0.5*exy, 0.5*exz],
        [0.5*exy, eyy, 0.5*eyz],
        [0.5*exz, 0.5*eyz, ezz],
    ])
    return E


def apply_strain(at: Atoms, eps_voigt: np.ndarray) -> Atoms:
    at2 = at.copy()
    E = strain_matrix_from_voigt(eps_voigt)
    F = np.eye(3) + E  # small strain approximation
    cell = at2.get_cell().array
    # ASE stores cell vectors as rows; apply right-acting F via cell * F^T
    new_cell = cell @ F.T
    at2.set_cell(new_cell, scale_atoms=True)
    return at2


def get_stress_GPa(at: Atoms, calc: CHGNetCalculator) -> np.ndarray:
    at.calc = calc
    s6 = at.get_stress(voigt=True)  # ASE sign convention: positive = compressive
    return np.array(s6, dtype=float) * GPA_PER_eV_per_A3


def relax_positions_only(at: Atoms, calc: CHGNetCalculator, tag: str,
                         fmax: float = 0.10, steps: int = 150) -> Atoms:
    at2 = at.copy()
    at2.calc = calc
    opt = BFGS(at2, logfile=f"strain_relax_{tag}.log")
    opt.run(fmax=fmax, steps=steps)
    return at2


def elastic_constants_via_stress_strain(at_eq: Atoms, calc: CHGNetCalculator,
                                        delta: float = 5e-3,
                                        relax_internal: bool = True,
                                        tag: str = "") -> np.ndarray:
    """Compute full 6x6 stiffness tensor C (GPa) by central differences.
    For each strain component j, compute stress at +delta and -delta after relaxing
    internal coordinates (positions-only) if relax_internal.
    """
    C = np.zeros((6, 6), dtype=float)
    for j in range(6):
        eps_p = np.zeros(6)
        eps_m = np.zeros(6)
        eps_p[j] = delta
        eps_m[j] = -delta

        at_p = apply_strain(at_eq, eps_p)
        at_m = apply_strain(at_eq, eps_m)

        if relax_internal:
            at_p = relax_positions_only(at_p, calc, tag=f"{tag}_j{j}_p")
            at_m = relax_positions_only(at_m, calc, tag=f"{tag}_j{j}_m")

        sig_p = get_stress_GPa(at_p, calc)
        sig_m = get_stress_GPa(at_m, calc)
        dsdE = (sig_p - sig_m) / (2.0 * delta)  # GPa per unit strain
        C[:, j] = dsdE
    return C


def vrh_from_C(C: np.ndarray):
    """Voigt–Reuss–Hill averages (B, G, E) from full 6x6 C (GPa)."""
    # Voigt averages
    C11, C22, C33 = C[0, 0], C[1, 1], C[2, 2]
    C12, C13, C23 = C[0, 1], C[0, 2], C[1, 2]
    C44, C55, C66 = C[3, 3], C[4, 4], C[5, 5]
    B_V = (C11 + C22 + C33 + 2.0 * (C12 + C13 + C23)) / 9.0
    G_V = (C11 + C22 + C33 - (C12 + C13 + C23) + 3.0 * (C44 + C55 + C66)) / 15.0

    # Reuss averages
    S = np.linalg.inv(C)
    S11, S22, S33 = S[0, 0], S[1, 1], S[2, 2]
    S12, S13, S23 = S[0, 1], S[0, 2], S[1, 2]
    S44, S55, S66 = S[3, 3], S[4, 4], S[5, 5]
    B_R = 1.0 / (S11 + S22 + S33 + 2.0 * (S12 + S13 + S23))
    G_R = 15.0 / (4.0 * (S11 + S22 + S33) - 4.0 * (S12 + S13 + S23) + 3.0 * (S44 + S55 + S66))

    B = 0.5 * (B_V + B_R)
    G = 0.5 * (G_V + G_R)
    E = 9.0 * B * G / (3.0 * B + G)
    nu = (3.0 * B - 2.0 * G) / (2.0 * (3.0 * B + G))
    return {
        "B_V": float(B_V), "G_V": float(G_V),
        "B_R": float(B_R), "G_R": float(G_R),
        "B": float(B), "G": float(G), "E": float(E), "nu": float(nu)
    }


@dataclass
class Result:
    label: str
    natoms: int
    composition: dict
    supercell: tuple
    a0_guess: float
    C_GPa: list
    vrh: dict
    E_GPa: float
    B_GPa: float
    G_GPa: float
    nu: float
    strain_delta: float
    notes: str


def main():
    seed_all(42)
    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    results = []

    # ----- Pure Al -----
    a0 = 4.05
    al = make_fcc_al(a0)
    al = supercell(al, (2, 2, 2))  # 32 atoms for speed
    ensure_min_distances(al)
    al_relaxed = relax_structure(al, calc, tag="Al_pure", variable_cell=True)
    write("POSCAR_Al_relaxed.vasp", al_relaxed, format="vasp")

    C_al = elastic_constants_via_stress_strain(al_relaxed, calc, delta=5e-3,
                                               relax_internal=True, tag="Al_pure")
    vrh_al = vrh_from_C(C_al)
    res_al = Result(
        label="Al",
        natoms=len(al_relaxed),
        composition=dict(Counter(al_relaxed.get_chemical_symbols())),
        supercell=(2, 2, 2),
        a0_guess=a0,
        C_GPa=C_al.tolist(),
        vrh=vrh_al,
        E_GPa=vrh_al["E"],
        B_GPa=vrh_al["B"],
        G_GPa=vrh_al["G"],
        nu=vrh_al["nu"],
        strain_delta=5e-3,
        notes="Positions relaxed for each strained state; cell fixed during elasticity."
    )
    results.append(asdict(res_al))

    # ----- Al–Cu ~1.85 at.% -----
    al2 = make_fcc_al(a0)
    al2 = supercell(al2, (3, 3, 3))  # 108 atoms
    # Substitute 2 Cu atoms -> 2/108 ~ 1.85 at.%
    al2 = substitute_cu(al2, n_cu=2, min_sep=7.5)
    ensure_min_distances(al2)
    al2_relaxed = relax_structure(al2, calc, tag="AlCu_1p85at", variable_cell=True)
    write("POSCAR_AlCu_relaxed.vasp", al2_relaxed, format="vasp")

    C_alcu = elastic_constants_via_stress_strain(al2_relaxed, calc, delta=5e-3,
                                                 relax_internal=True, tag="AlCu_1p85at")
    vrh_alcu = vrh_from_C(C_alcu)
    res_alcu = Result(
        label="Al-1.85at%Cu",
        natoms=len(al2_relaxed),
        composition=dict(Counter(al2_relaxed.get_chemical_symbols())),
        supercell=(3, 3, 3),
        a0_guess=a0,
        C_GPa=C_alcu.tolist(),
        vrh=vrh_alcu,
        E_GPa=vrh_alcu["E"],
        B_GPa=vrh_alcu["B"],
        G_GPa=vrh_alcu["G"],
        nu=vrh_alcu["nu"],
        strain_delta=5e-3,
        notes="Positions relaxed for each strained state; cell fixed during elasticity. Cu sites far-separated."
    )
    results.append(asdict(res_alcu))

    # Summaries
    # Sort to have Al first
    results_sorted = sorted(results, key=lambda x: x["label"])
    with open("results.json", "w") as f:
        json.dump({"results": results_sorted}, f, indent=2)

    # CSV summary
    import csv
    with open("summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "natoms", "E_GPa", "B_GPa", "G_GPa", "nu"])
        for r in results_sorted:
            w.writerow([r["label"], r["natoms"], f"{r['E_GPa']:.3f}", f"{r['B_GPa']:.3f}", f"{r['G_GPa']:.3f}", f"{r['nu']:.4f}"])

    # Human summary
    # Extract E values
    E_al = next(r["E_GPa"] for r in results_sorted if r["label"] == "Al")
    E_alcu = next(r["E_GPa"] for r in results_sorted if r["label"].startswith("Al-"))
    pct = 100.0 * (E_alcu - E_al) / E_al

    print("Elastic properties (Voigt–Reuss–Hill):")
    for r in results_sorted:
        print(f"- {r['label']}: E = {r['E_GPa']:.2f} GPa, B = {r['B_GPa']:.2f} GPa, G = {r['G_GPa']:.2f} GPa, nu = {r['nu']:.3f}")
    print(f"\nDelta E: {E_alcu - E_al:.2f} GPa ({pct:.1f}%)")
    print("Conclusion: Cu increases E by at least 5%? ", "YES" if pct >= 5.0 else "NO")


if __name__ == "__main__":
    main()
