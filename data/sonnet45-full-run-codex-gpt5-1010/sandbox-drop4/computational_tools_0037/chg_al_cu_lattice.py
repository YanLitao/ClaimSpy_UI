import os
import json
import csv
import numpy as np

# Reproducibility
np.random.seed(42)

from ase.build import bulk
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator


def make_al_supercell(a0=4.05, rep=(3, 3, 3)):
    """Build a conventional fcc Al cell and make a supercell.

    - a0: initial lattice parameter [Å]
    - rep: repetition of the conventional cell; (3,3,3) gives 4*27 = 108 atoms
    """
    al_conv = bulk("Al", "fcc", a=a0, cubic=True)
    sc = al_conv.repeat(rep)
    return sc


def substitute_cu(atoms, n_cu, seed=42):
    """Randomly substitute n_cu Al atoms with Cu. Returns new Atoms and indices.
    """
    if n_cu <= 0:
        return atoms.copy(), []
    rng = np.random.default_rng(seed)
    idxs = rng.choice(len(atoms), size=n_cu, replace=False)
    out = atoms.copy()
    symbols = out.get_chemical_symbols()
    for i in idxs:
        symbols[i] = "Cu"
    out.set_chemical_symbols(symbols)
    return out, idxs.tolist()


def min_interatomic_distance(atoms):
    # Simple O(N^2) check; acceptable for ~100 atoms
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()
    nat = len(atoms)
    dmin = float("inf")
    for i in range(nat - 1):
        for j in range(i + 1, nat):
            rij = pos[j] - pos[i]
            if any(pbc):
                # Minimum image convention
                rij = rij - np.round(np.linalg.solve(cell.T, rij)) @ cell
            d = np.linalg.norm(rij)
            if d < dmin:
                dmin = d
    return dmin


def relax_structure(atoms, calc, label, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300):
    """Two-stage relaxation: positions-only, then hydrostatic cell + positions.
    Returns relaxed atoms and a metrics dict.
    """
    atoms = atoms.copy()
    atoms.set_calculator(calc)

    # Guard against pathological starting states
    e0 = atoms.get_potential_energy()
    f0 = atoms.get_forces()
    if not np.isfinite(e0) or not np.isfinite(f0).all():
        raise RuntimeError("Non-finite energy/forces before optimization; check initial structure.")

    # Stage 1: positions-only (keep cell fixed)
    opt1 = BFGS(atoms, logfile=f"relax1_{label}.log")
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)

    e1 = atoms.get_potential_energy()
    f1 = atoms.get_forces()
    if not np.isfinite(e1) or not np.isfinite(f1).all():
        raise RuntimeError("Non-finite energy/forces after stage 1; aborting.")

    # Stage 2: allow isotropic cell strains (hydrostatic) + positions
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile=f"relax2_{label}.log")
    opt2.run(fmax=fmax_fine, steps=max_steps)

    ef = atoms.get_potential_energy()
    ff = atoms.get_forces()
    if not np.isfinite(ef) or not np.isfinite(ff).all():
        raise RuntimeError("Non-finite energy/forces after stage 2; aborting.")

    metrics = {
        "E_total_eV": float(ef),
        "E_per_atom_eV": float(ef / len(atoms)),
    }
    return atoms, metrics


def equivalent_cubic_a_from_volume(atoms, n_conv_cells):
    """Compute equivalent cubic lattice parameter a from the supercell volume.
    a = (V_supercell / n_conv_cells)^(1/3)
    """
    V = atoms.get_volume()
    return float((V / n_conv_cells) ** (1.0 / 3.0))


def main():
    # Fix seeds
    try:
        import torch

        torch.manual_seed(42)
    except Exception:
        pass

    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Build a 3x3x3 supercell of the conventional fcc cell: 4*27 = 108 atoms
    rep = (3, 3, 3)
    n_conv_cells = rep[0] * rep[1] * rep[2]
    base = make_al_supercell(a0=4.05, rep=rep)

    # Sanity: minimum distance should be reasonable
    dmin0 = min_interatomic_distance(base)
    if dmin0 < 0.7:
        # Expand slightly if truly pathological (very unlikely here)
        base.set_cell(base.get_cell() * (0.7 / dmin0), scale_atoms=True)

    cases = [(0, "Al")] + [(n, f"AlCu_{n}Cu") for n in range(1, 6)]

    results = []

    for n_cu, label in cases:
        if n_cu == 0:
            atoms = base.copy()
            subs = []
        else:
            # Different seeds to diversify patterns but keep reproducible
            atoms, subs = substitute_cu(base, n_cu, seed=42 + n_cu)

        # Relax
        relaxed, metrics = relax_structure(atoms, calc, label, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300)

        # Lattice parameter from volume per conventional cell
        a_eq = equivalent_cubic_a_from_volume(relaxed, n_conv_cells=n_conv_cells)

        # Save POSCAR
        write(f"POSCAR_{label}.vasp", relaxed, vasp5=True, direct=True, sort=True)

        # Compose record
        n_cu_final = sum(1 for s in relaxed.get_chemical_symbols() if s == "Cu")
        x_cu_pct = 100.0 * n_cu_final / len(relaxed)
        rec = {
            "label": label,
            "natoms": len(relaxed),
            "n_Cu": int(n_cu_final),
            "x_Cu_pct": float(x_cu_pct),
            "a_eq_A": float(a_eq),
            "E_per_atom_eV": float(metrics["E_per_atom_eV"]),
        }
        results.append(rec)

    # Reference lattice parameter (pure Al)
    a0 = next(r["a_eq_A"] for r in results if r["n_Cu"] == 0)
    for r in results:
        r["delta_a_over_a0_pct"] = 100.0 * (r["a_eq_A"] / a0 - 1.0)

    # Write machine-readable outputs
    with open("results.json", "w") as f:
        json.dump(results, f, indent=2)

    with open("summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "natoms", "n_Cu", "x_Cu_pct", "a_eq_A", "delta_a_over_a0_pct", "E_per_atom_eV"])
        for r in results:
            w.writerow([
                r["label"],
                r["natoms"],
                r["n_Cu"],
                f"{r['x_Cu_pct']:.3f}",
                f"{r['a_eq_A']:.6f}",
                f"{r['delta_a_over_a0_pct']:.4f}",
                f"{r['E_per_atom_eV']:.6f}",
            ])

    # Human summary
    print("Computed equilibrium lattice parameters (equivalent cubic a):")
    print("label        natoms  x_Cu[%]   a_eq[Å]    Δa/a0[%]")
    for r in results:
        print(
            f"{r['label']:<12} {r['natoms']:>6}  {r['x_Cu_pct']:>6.3f}  "
            f"{r['a_eq_A']:>8.5f}  {r['delta_a_over_a0_pct']:>8.4f}"
        )

    # Explicitly report the requested items
    pure = next(r for r in results if r["n_Cu"] == 0)
    print("\n(1) Pure Al lattice parameter a0 [Å]:", f"{pure['a_eq_A']:.6f}")

    print("(2) Al-Cu lattice parameters by concentration:")
    for r in results:
        if r["n_Cu"] > 0:
            print(
                f"    {r['label']}: x_Cu={r['x_Cu_pct']:.3f}% -> a={r['a_eq_A']:.6f} Å"
            )

    print("(3) Relative change Δa/a0 [%]:")
    for r in results:
        if r["n_Cu"] > 0:
            print(
                f"    {r['label']}: x_Cu={r['x_Cu_pct']:.3f}% -> Δa/a0={r['delta_a_over_a0_pct']:.4f}%"
            )


if __name__ == "__main__":
    # Ensure all outputs land in CWD and not elsewhere
    os.chdir(os.getcwd())
    main()

