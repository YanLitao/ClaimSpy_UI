#!/usr/bin/env python3
"""
Compute the bulk modulus of elemental Ag (fcc) using CHGNet via an
energy–volume (E–V) curve and a Birch–Murnaghan EOS fit.

Workflow (single script, single CHGNet load):
- Build fcc Ag near experimental lattice (a≈4.08 Å)
- Two-stage relaxation (positions → hydrostatic cell)
- Sample E(V) by isotropic cell scaling around V0 (single-point energies)
- Fit Birch–Murnaghan EOS to extract B0 (reported in GPa)
- Write: results.json, summary.csv, e_v_curve.csv, Ag_relaxed.vasp, eos_fit.png

Notes:
- Uses only current working directory for intermediates.
- Loads .env (for MP_API_KEY) but does not require MP queries.
- Avoids interactive calls; plotting is saved to file.
"""

import os
import json
import csv
import math
import random
import datetime as dt
from pathlib import Path

import numpy as np
from dotenv import load_dotenv

from ase.build import bulk
from ase.io import write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.eos import EquationOfState

from chgnet.model import CHGNet, CHGNetCalculator


# Constants
EV_PER_A3_TO_GPA = 160.21766208  # 1 eV/Å^3 = 160.21766208 GPa


def safe_min_distance(atoms) -> float:
    """Return the minimum interatomic distance with PBC MIC, excluding zeros."""
    dmat = atoms.get_all_distances(mic=True)
    # Mask zero diagonals
    dmat[dmat == 0] = np.inf
    return float(np.min(dmat))


def ensure_reasonable_geometry(atoms, min_allowed=0.7):
    """If any pair is closer than min_allowed Å, uniformly expand the cell until safe."""
    dmin = safe_min_distance(atoms)
    if np.isfinite(dmin) and dmin < min_allowed:
        # scale factor to bring dmin to min_allowed with small safety margin
        factor = (min_allowed / max(dmin, 1e-6)) * 1.05
        # isotropic scaling of cell and positions
        cell = atoms.get_cell()
        atoms.set_cell(cell * factor, scale_atoms=True)


def assert_finite_energy_forces(atoms):
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not np.isfinite(e) or not np.isfinite(f).all():
        raise RuntimeError("Non-finite energy/forces encountered.")
    return e


def relax_structure(atoms, calc, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300):
    atoms = atoms.copy()
    atoms.set_calculator(calc)

    ensure_reasonable_geometry(atoms)

    # Stage 1: positions-only
    opt1 = BFGS(atoms, logfile="relax1.log")
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)
    assert_finite_energy_forces(atoms)

    # Stage 2: hydrostatic cell relax
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile="relax2.log")
    opt2.run(fmax=fmax_fine, steps=max_steps)
    assert_finite_energy_forces(atoms)

    return atoms


def compute_e_v_curve(atoms_eq, calc, scales):
    volumes = []
    energies = []
    for s in scales:
        at = atoms_eq.copy()
        factor = s ** (1.0 / 3.0)  # convert volume scale to length scale
        at.set_cell(at.get_cell() * factor, scale_atoms=True)
        at.set_calculator(calc)
        e = at.get_potential_energy()
        v = at.get_volume()
        volumes.append(v)
        energies.append(e)
    return np.array(volumes), np.array(energies)


def fit_birch_murnaghan(volumes, energies):
    eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
    v0, e0, B = eos.fit()  # B in eV/Å^3
    B_gpa = float(B * EV_PER_A3_TO_GPA)
    return float(v0), float(e0), float(B), B_gpa, eos


def write_csv(path, rows, header):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        if header:
            w.writerow(header)
        w.writerows(rows)


def main():
    # Environment and reproducibility
    load_dotenv()  # For MP_API_KEY if needed in future
    np.random.seed(42)
    random.seed(42)

    # Build initial fcc Ag (conventional cubic cell)
    a0 = 4.08  # Å, room-temperature experimental starting guess
    atoms = bulk('Ag', 'fcc', a=a0, cubic=True)

    # Load CHGNet once and reuse
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Relax structure (positions then hydrostatic cell)
    atoms_eq = relax_structure(atoms, calc, fmax_coarse=0.2, fmax_fine=0.05, max_steps=300)

    # Persist the relaxed structure
    write('Ag_relaxed.vasp', atoms_eq, format='vasp', direct=True, vasp5=True)

    # Basic relaxed properties
    e_eq = atoms_eq.get_potential_energy()
    v_eq = atoms_eq.get_volume()
    natoms = len(atoms_eq)
    lengths = atoms_eq.get_cell().lengths()
    angles = atoms_eq.get_cell().angles()

    # Energy–volume sampling: symmetric about V0 with moderate spread
    # Use single-point energies (no internal DOFs for fcc Ag)
    scales = np.array([0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06])
    volumes, energies = compute_e_v_curve(atoms_eq, calc, scales)

    # Fit EOS
    v0, e0, B_eVa3, B_gpa, eos = fit_birch_murnaghan(volumes, energies)

    # Plot EOS fit
    try:
        eos.plot('eos_fit.png', show=False)
    except Exception:
        # Matplotlib may be unavailable; silently skip plotting
        pass

    # Write datasets
    # E–V points
    e_v_rows = [(float(s), float(v), float(e)) for s, v, e in zip(scales, volumes, energies)]
    write_csv('e_v_curve.csv', e_v_rows, header=['scale(V/V0_relaxed)', 'volume(Å^3)', 'energy(eV)'])

    # Summary CSV (one-liners)
    summary_rows = [
        ['material', 'Ag (fcc)'],
        ['n_atoms_cell', natoms],
        ['relaxed_a(Å)', f"{lengths[0]:.6f}"],
        ['relaxed_b(Å)', f"{lengths[1]:.6f}"],
        ['relaxed_c(Å)', f"{lengths[2]:.6f}"],
        ['relaxed_alpha(deg)', f"{angles[0]:.6f}"],
        ['relaxed_beta(deg)', f"{angles[1]:.6f}"],
        ['relaxed_gamma(deg)', f"{angles[2]:.6f}"],
        ['relaxed_volume_cell(Å^3)', f"{v_eq:.8f}"],
        ['relaxed_energy_cell(eV)', f"{e_eq:.8f}"],
        ['EOS', 'Birch-Murnaghan 3rd order'],
        ['V0_cell(Å^3)', f"{v0:.8f}"],
        ['E0_cell(eV)', f"{e0:.8f}"],
        ['B0(eV/Å^3)', f"{B_eVa3:.8f}"],
        ['B0(GPa)', f"{B_gpa:.6f}"],
    ]
    write_csv('summary.csv', summary_rows, header=['property', 'value'])

    # Results JSON
    results = {
        'material': 'Ag',
        'structure': 'fcc',
        'n_atoms_cell': int(natoms),
        'relaxed': {
            'a_b_c_A': [float(x) for x in lengths],
            'alpha_beta_gamma_deg': [float(x) for x in angles],
            'volume_cell_A3': float(v_eq),
            'energy_cell_eV': float(e_eq),
        },
        'eos': {
            'type': 'birchmurnaghan',
            'volumes_A3': [float(v) for v in volumes],
            'energies_eV': [float(e) for e in energies],
            'scales_V_over_Veq': [float(s) for s in scales],
            'V0_cell_A3': float(v0),
            'E0_cell_eV': float(e0),
            'B0_eV_per_A3': float(B_eVa3),
            'B0_GPa': float(B_gpa),
        },
        'units': {
            'length': 'Å',
            'energy': 'eV',
            'volume': 'Å^3',
            'bulk_modulus': 'GPa',
        },
        'method': 'CHGNet (MLIP) + ASE EOS fit',
        'timestamp_utc': dt.datetime.utcnow().isoformat() + 'Z',
    }
    with open('results.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Human summary
    v0_pa = v0 / natoms
    v_eq_pa = v_eq / natoms
    print("=== CHGNet E–V EOS Fit for Ag (fcc) ===")
    print(f"Relaxed a (Å): {lengths[0]:.6f}; cell volume (Å^3): {v_eq:.6f} ({v_eq_pa:.6f} Å^3/atom)")
    print(f"EOS V0 (Å^3): {v0:.6f} ({v0_pa:.6f} Å^3/atom)")
    print(f"Bulk modulus B0: {B_gpa:.2f} GPa")
    print("Files written: Ag_relaxed.vasp, e_v_curve.csv, summary.csv, results.json, eos_fit.png")


if __name__ == '__main__':
    main()

