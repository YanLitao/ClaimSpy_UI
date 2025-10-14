import os
import json
import csv
import math
import numpy as np

# Seed control
np.random.seed(42)
try:
    import torch
    torch.manual_seed(42)
except Exception:
    pass

# Load MP key from .env if present (not strictly needed here)
try:
    from dotenv import load_dotenv
    load_dotenv(dotenv_path=".env")
except Exception:
    pass

from ase.build import bulk
from ase.io import write
from ase.optimize import BFGS
from ase.eos import EquationOfState
from ase.atoms import Atoms

from chgnet.model import CHGNet, CHGNetCalculator

EV_PER_A3_TO_GPA = 160.21766208


def build_ag_supercell(a0=4.09, rep=(2, 2, 2)) -> Atoms:
    # fcc, primitive cell (1 atom), then replicate to get 8 atoms total
    ag_prim = bulk('Ag', 'fcc', a=a0)
    ag = ag_prim.repeat(rep)
    ag.pbc = True
    return ag


def build_agcu_supercell(a0=4.09, rep=(2, 2, 2), x_cu=0.125) -> Atoms:
    ag = build_ag_supercell(a0=a0, rep=rep)
    n = len(ag)
    n_cu = max(1, int(round(x_cu * n)))
    # Guard to hit exactly 12.5% when n=8
    if n == 8 and abs(x_cu - 0.125) < 1e-6:
        n_cu = 1
    # Replace first n_cu atoms with Cu; with 8 atoms and 12.5%, this is 1 Cu
    symbols = ag.get_chemical_symbols()
    for i in range(n_cu):
        symbols[i] = 'Cu'
    ag.set_chemical_symbols(symbols)
    return ag


def positions_relax(atoms: Atoms, calc: CHGNetCalculator, fmax=0.1, steps=150):
    atoms.set_calculator(calc)
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)


def energy_vs_volume(atoms: Atoms, calc: CHGNetCalculator, lin_scales=None, fmax=0.1):
    if lin_scales is None:
        lin_scales = np.linspace(0.96, 1.04, 9)
    base_cell = atoms.cell.array.copy()
    base_pos = atoms.get_positions().copy()

    volumes = []
    energies = []
    structures = []  # store the relaxed structures at sampled volumes

    for k in lin_scales:
        a = atoms.copy()
        # Isotropic linear scaling of the cell; scale positions with cell
        a.set_cell(base_cell * k, scale_atoms=True)
        a.set_positions(a.get_positions())
        a.pbc = True
        a.set_calculator(calc)
        # Relax positions only at fixed cell
        opt = BFGS(a, logfile=None)
        opt.run(fmax=fmax, steps=150)
        e = a.get_potential_energy()
        v = a.get_volume()
        if not (np.isfinite(e) and np.isfinite(v)):
            raise RuntimeError("Non-finite energy/volume encountered during EOS sampling.")
        volumes.append(v)
        energies.append(e)
        structures.append(a)

    return np.array(volumes), np.array(energies), structures


def fit_eos(volumes, energies, label, plot_filename=None):
    eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
    v0, e0, b = eos.fit()  # b in eV/Å^3
    if plot_filename is not None:
        try:
            eos.plot(plot_filename)
        except Exception:
            pass
    B_gpa = float(b) * EV_PER_A3_TO_GPA
    return float(v0), float(e0), float(b), B_gpa


def main():
    # Load CHGNet model once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Build structures
    ag = build_ag_supercell(a0=4.09, rep=(2, 2, 2))
    agcu = build_agcu_supercell(a0=4.09, rep=(2, 2, 2), x_cu=0.125)

    # Quick pre-check energies/forces to guard against NaNs
    for name, a in [("Ag", ag), ("Ag0.875Cu0.125", agcu)]:
        a.set_calculator(calc)
        e = a.get_potential_energy()
        f = a.get_forces()
        if not (np.isfinite(e) and np.isfinite(f).all()):
            raise RuntimeError(f"Non-finite energy/forces for {name} before optimization.")

    # EOS sampling with positions-only relax across volumes
    scales = np.linspace(0.96, 1.04, 9)

    v_ag, e_ag, structs_ag = energy_vs_volume(ag, calc, lin_scales=scales, fmax=0.1)
    v0_ag, e0_ag, b_ag_eva3, b_ag_gpa = fit_eos(v_ag, e_ag, 'Ag', plot_filename='eos_ag.png')

    v_agcu, e_agcu, structs_agcu = energy_vs_volume(agcu, calc, lin_scales=scales, fmax=0.1)
    v0_agcu, e0_agcu, b_agcu_eva3, b_agcu_gpa = fit_eos(v_agcu, e_agcu, 'Ag0.875Cu0.125', plot_filename='eos_agcu.png')

    # Pick the sampled structure nearest to v0 for writing POSCARs
    idx_ag = int(np.abs(v_ag - v0_ag).argmin())
    idx_agcu = int(np.abs(v_agcu - v0_agcu).argmin())
    write('POSCAR_Ag', structs_ag[idx_ag], format='vasp')
    write('POSCAR_Ag0.875Cu0.125', structs_agcu[idx_agcu], format='vasp')

    # Summaries
    results = {
        'units': {
            'bulk_modulus': 'GPa',
            'volume': 'Å^3',
            'energy': 'eV',
            'b_eos': 'eV/Å^3'
        },
        'structures': [
            {
                'label': 'Ag',
                'composition': 'Ag',
                'n_atoms': len(ag),
                'v0': v0_ag,
                'e0': e0_ag,
                'b_eva3': b_ag_eva3,
                'bulk_GPa': b_ag_gpa,
            },
            {
                'label': 'Ag0.875Cu0.125',
                'composition': 'Ag7Cu (8-atom fcc supercell with 1 Cu; 12.5 at.% Cu)',
                'n_atoms': len(agcu),
                'v0': v0_agcu,
                'e0': e0_agcu,
                'b_eva3': b_agcu_eva3,
                'bulk_GPa': b_agcu_gpa,
            }
        ]
    }

    # Percent change
    pct_change = (b_agcu_gpa - b_ag_gpa) / b_ag_gpa * 100.0

    # Write JSON and CSV summaries
    with open('results.json', 'w') as f:
        json.dump({**results, 'percent_change_bulk_modulus': pct_change}, f, indent=2)

    with open('summary.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Label', 'Bulk_modulus_GPa', 'v0_A3', 'e0_eV', 'b_eVA3'])
        w.writerow(['Ag', f"{b_ag_gpa:.3f}", f"{v0_ag:.3f}", f"{e0_ag:.6f}", f"{b_ag_eva3:.6f}"])
        w.writerow(['Ag0.875Cu0.125', f"{b_agcu_gpa:.3f}", f"{v0_agcu:.3f}", f"{e0_agcu:.6f}", f"{b_agcu_eva3:.6f}"])
        w.writerow(['Percent change (alloy vs Ag)', f"{pct_change:.2f}", '', '', ''])

    # Print concise human summary for CLI
    print("FoM: Bulk modulus B (GPa)")
    print(f"Ag: {b_ag_gpa:.2f} GPa")
    print(f"Ag0.875Cu0.125: {b_agcu_gpa:.2f} GPa")
    print(f"Percent change: {pct_change:.2f}% (alloy vs Ag)")


if __name__ == '__main__':
    main()
