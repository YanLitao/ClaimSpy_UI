#!/usr/bin/env python3
"""
CALPHAD equilibrium for Cu–Al–Li at 298 K using available TDB(s).

Workflow:
- Load database from the TnETDBDB directory
- Verify required elements are present (AL, CU, LI, VA)
- Build phase list by introspecting the database
- Compute equilibrium at 298 K for X(CU)=0.30, X(AL)=0.30 (X(LI)=0.40 implied)
- Derive phase fractions; list stable phases at T=298 K
- Write results.json and results.csv, and print a concise summary

All outputs are written in the current working directory.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pycalphad import Database, equilibrium, variables as v


def main():
    # Fixed inputs per the prompt
    T = 298.0  # K
    P = 101325.0  # Pa
    x_cu = 0.30
    x_al = 0.30
    # X(LI) implied by N−1 rule for mole fractions in substitutional solution

    # Locate TDB directory (use the first available per policy)
    tdb_dirs = [
        Path('/Users/delip/play/k2code/TnETDBDB'),
        Path('/Users/delip/play/miniclaimspy-service/TnETDBDB'),
        Path('/Users/delip/play/miniclaimspy/TnETDBDB'),
    ]
    tdb_path = None
    for d in tdb_dirs:
        if d.is_dir():
            # Prefer any .tdb/.TDB file; pick the first deterministically
            cands = sorted(list(d.glob('*.tdb')) + list(d.glob('*.TDB')))
            if cands:
                tdb_path = cands[0]
                break

    out_json = Path('results.json')
    out_csv = Path('results.csv')

    def write_error(reason: str):
        payload = {
            'status': 'error',
            'reason': reason,
            'tdb': str(tdb_path) if tdb_path else None,
        }
        out_json.write_text(json.dumps(payload, indent=2))
        # Empty CSV with headers for consistency
        pd.DataFrame(columns=['T', 'phase', 'fraction']).to_csv(out_csv, index=False)
        print(f"ERROR: {reason}")
        return 2

    if tdb_path is None:
        return write_error('No TDB file found in any TnETDBDB directory.')

    try:
        dbf = Database(str(tdb_path))
    except Exception as exc:
        return write_error(f'Failed to load TDB: {exc}')

    required = {'AL', 'CU', 'LI', 'VA'}
    # pycalphad Database exposes elements in dbf.elements
    db_elems = set(getattr(dbf, 'elements', []))
    missing = required - db_elems
    if missing:
        return write_error(
            f"Database '{tdb_path.name}' missing required elements: {sorted(missing)}; "
            f"available: {sorted(db_elems)}"
        )

    # Components and phase selection
    comps = ['AL', 'CU', 'LI', 'VA']

    # Use all phases defined in the database (generic approach per guidance)
    phases = sorted(dbf.phases.keys())
    if not phases:
        return write_error('No phases selected; introspection returned an empty phase list.')

    # Conditions: enforce N−1 mole fraction constraints for substitutional species
    conds = {
        v.T: T,
        v.P: P,
        v.N: 1.0,
        v.X('CU'): x_cu,
        v.X('AL'): x_al,
    }
    # Composition constraints follow the N−1 rule for a ternary (two X() constraints)

    try:
        eq = equilibrium(dbf, comps, phases, conds, verbose=False)
    except Exception as exc:
        return write_error(f'equilibrium() failed: {exc}')

    # Phase fractions from NP; sum vertices if present
    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')
    if 'phase' in pf.dims:
        tot = pf.sum('phase')
        # guard against divide-by-zero
        with np.errstate(divide='ignore', invalid='ignore'):
            phase_frac = pf / tot

        # Extract at the single condition (T=298 K); flatten to a dict
        # phase dimension remains; other dims (P, comps, etc.) are size-1
        pf_298 = phase_frac.squeeze()
        phases_list = list(pf_298['phase'].values)
        frac_vals = np.array(pf_298.values, dtype=float)

        # Determine stable phases with a practical threshold
        stable_mask = frac_vals > 1e-3
        stable_phases = [ph for ph, m in zip(phases_list, stable_mask) if m]
        stable_fracs = [float(f) for f, m in zip(frac_vals, stable_mask) if m]

        # Prepare outputs
        records = []
        for ph, fr in zip(phases_list, frac_vals):
            records.append({'T': T, 'phase': ph, 'fraction': float(fr)})
        df = pd.DataFrame.from_records(records)
        df.to_csv(out_csv, index=False)

        summary = {
            'status': 'ok',
            'tdb': str(tdb_path),
            'T_K': T,
            'P_Pa': P,
            'composition_atfrac': {'CU': x_cu, 'AL': x_al, 'LI': 1.0 - x_cu - x_al},
            'stable_phases': [
                {'phase': ph, 'fraction': fr} for ph, fr in sorted(
                    zip(stable_phases, stable_fracs), key=lambda x: -x[1]
                )
            ],
            'single_phase': (sum(fr > 1e-3 for fr in stable_fracs) == 1),
        }
    else:
        # No 'phase' dimension returned; cannot infer phase fractions
        df = pd.DataFrame(columns=['T', 'phase', 'fraction'])
        df.to_csv(out_csv, index=False)
        summary = {
            'status': 'ok',
            'warning': "Equilibrium result has no 'phase' dimension; this TDB likely lacks parameterization for Cu–Al–Li.",
            'tdb': str(tdb_path),
            'T_K': T,
            'P_Pa': P,
            'composition_atfrac': {'CU': x_cu, 'AL': x_al, 'LI': 1.0 - x_cu - x_al},
            'stable_phases': [],
            'single_phase': False,
        }
    out_json.write_text(json.dumps(summary, indent=2))

    # Human summary
    print(f"TDB: {tdb_path.name}")
    print(f"Conditions: T={T:.1f} K, P={P:.0f} Pa, X(CU)={x_cu:.2f}, X(AL)={x_al:.2f}, X(LI)={1.0 - x_cu - x_al:.2f}")
    if summary.get('warning'):
        print(f"Warning: {summary['warning']}")
        print('No phase fractions available from this calculation with the given database.')
    elif not summary['stable_phases']:
        print('At 298 K, no stable phases identified (zero total phase amount).')
    else:
        print("Stable phases at 298 K (>1e-3 fraction):")
        for item in summary['stable_phases']:
            print(f"  - {item['phase']}: {item['fraction']:.6f}")
        if summary['single_phase']:
            print('Conclusion: Single-phase equilibrium at this composition and T.')
        else:
            print('Conclusion: Multi-phase equilibrium (decomposition into the above phases).')

    return 0


if __name__ == '__main__':
    sys.exit(main())
