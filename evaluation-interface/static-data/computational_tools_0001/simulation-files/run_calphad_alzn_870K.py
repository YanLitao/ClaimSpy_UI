#!/usr/bin/env python3
"""
CALPHAD equilibrium calculation for Al–Zn at 870 K, X(AL)=0.20 (at.%).

Workflow:
- Load COST507 database from TnETDBDB
- Select relevant phases from DB (containing AL/ZN)
- Solve equilibrium at T=870 K, P=1 atm, N=1 mol, X(AL)=0.20
- Compute phase fractions from NP, threshold for stability at 1e-3
- Emit results.json and results.csv
- Print concise human-readable summary

Notes:
- Follows the one-script, single-run guidance.
- Uses pycalphad to introspect the database; no grepping TDB files.
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from dotenv import load_dotenv  # for MP API if needed; not used here
except Exception:
    load_dotenv = None

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf: Database, phase: str, comps):
    const = dbf.phases[phase].constituents
    flat = {c for subl in const for c in subl}
    # Map Species objects to their element names (e.g., 'AL', 'ZN')
    elems = {getattr(c, 'name', str(c)) for c in flat}
    return bool(set(comps) & elems)


def main():
    # Load .env if present (e.g., for MP_API_KEY); no external queries here
    if load_dotenv is not None:
        env_path = Path('.env')
        if env_path.exists():
            load_dotenv(env_path)

    # Database path (COST507, under TnETDBDB). Read-only.
    db_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    if not db_path.exists():
        raise FileNotFoundError(f"TDB database not found: {db_path}")

    dbf = Database(str(db_path))

    # Components: include VA if present in DB
    target_comps = ['AL', 'ZN', 'VA']
    comps = [c for c in target_comps if c in dbf.elements]
    if ('AL' not in comps) or ('ZN' not in comps):
        raise ValueError(f"Database does not include required elements AL/ZN; has {sorted(dbf.elements)}")

    # Phases containing AL or ZN; sorted for reproducibility
    phases = sorted([p for p in dbf.phases.keys() if phase_has_any(dbf, p, ['AL', 'ZN'])])
    assert phases, "No phases selected; did you build the phases list?"

    # Conditions: N-1 mole fraction constraints for a binary system
    # Temperature in Kelvin, pressure in Pa
    T = 870.0
    P = 101325.0
    X_AL = 0.20  # 20 at.% Al => 80 at.% Zn
    conds = {v.T: T, v.P: P, v.N: 1.0, v.X('AL'): X_AL}

    # (Guard omitted; composition constraints are set explicitly above.)

    # Solve equilibrium
    eq = equilibrium(dbf, comps, phases, conds, output='NP')

    # Compute phase fractions from vertex labeling in this pycalphad version.
    # For each vertex, look up phase name and moles NP; accumulate per phase.
    np_da = eq['NP']  # dims: (N, P, T, X_AL, vertex)
    ph_da = eq['Phase']  # dims: (N, P, T, X_AL, vertex)

    # Extract along the single condition point (N,P,T,X_AL all len=1)
    np_vals = np_da.values.reshape(-1)
    ph_vals = ph_da.values.reshape(-1)

    phase_moles = {}
    for ph_label, np_mol in zip(ph_vals, np_vals):
        ph_name = str(ph_label).strip()
        if not ph_name:
            continue
        phase_moles[ph_name] = phase_moles.get(ph_name, 0.0) + float(np_mol)

    tot_moles = sum(phase_moles.values())
    if tot_moles <= 0:
        raise RuntimeError("Total phase moles computed as zero or negative; check equilibrium result.")

    results = [{'phase': ph, 'fraction': m / tot_moles} for ph, m in phase_moles.items()]

    # Stability threshold per guidance
    stable = [r for r in results if r['fraction'] > 1e-3]
    stable_sorted = sorted(stable, key=lambda r: r['fraction'], reverse=True)

    # Identify liquid if present
    liquid_frac = 0.0
    for r in results:
        if r['phase'].upper().startswith('LIQUID'):
            liquid_frac = r['fraction']
            break

    # Write outputs
    out_json = Path('alzn_870K_results.json')
    out_csv = Path('alzn_870K_results.csv')

    payload = {
        'database': str(db_path),
        'components': comps,
        'phases_considered': phases,
        'conditions': {'T_K': T, 'P_Pa': P, 'X_AL': X_AL, 'N_mol': 1.0},
        'phase_fractions': results,
        'stable_phases_threshold_1e-3': stable_sorted,
        'liquid_present': liquid_frac > 1e-6,
        'liquid_fraction': liquid_frac,
    }
    with out_json.open('w') as f:
        json.dump(payload, f, indent=2)

    # CSV of all phases; write stable flag
    df = pd.DataFrame(results)
    df['stable_gt_1e-3'] = df['fraction'] > 1e-3
    df.to_csv(out_csv, index=False)

    # Human-readable summary
    print("Al–Zn CALPHAD @ 870 K, X(AL)=0.20 (P=1 atm)")
    print(f"Database: {db_path}")
    if stable_sorted:
        print("Stable phases (fraction, >1e-3):")
        for r in stable_sorted:
            print(f"  - {r['phase']}: {r['fraction']:.6f}")
    else:
        print("No phases exceed the 1e-3 fraction threshold; showing top contributors:")
        top = sorted(results, key=lambda r: r['fraction'], reverse=True)[:5]
        for r in top:
            print(f"  - {r['phase']}: {r['fraction']:.6f}")

    state = "Liquid present" if (liquid_frac > 1e-6) else "All solid"
    print(f"State assessment: {state} (LIQUID fraction = {liquid_frac:.6f})")


if __name__ == '__main__':
    main()
