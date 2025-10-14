#!/usr/bin/env python3
"""
Compute the equilibrium phase constitution for Cu80Li10Al10 (at.%) at 298 K and 1 atm
using a local TDB from the TnETDBDB directory, following the CALPHAD.md workflow.

Outputs:
- results.json: machine-readable summary
- results.csv: phases and fractions
Prints a concise human summary to stdout.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import xarray as xr
from pycalphad import Database, equilibrium, variables as v


def find_local_tdbs(tdb_dir: Path):
    tdbs = []
    if not tdb_dir.exists():
        return tdbs
    for p in sorted(tdb_dir.iterdir()):
        if p.suffix.lower() == ".tdb":
            tdbs.append(p)
    return tdbs


def db_covers_elements(db: Database, needed):
    # pycalphad Database has .elements including VA; ensure needed subset present
    try:
        present = set(db.elements)
    except Exception:
        # Fallback: infer from species strings
        present = set()
        for sp in getattr(db, 'species', []):
            sym = str(sp).strip().upper()
            # crude parse: species like AL, CU, LI, AL+, AL-; take letters
            if sym.isalpha():
                present.add(sym)
    return set(map(str.upper, needed)).issubset({e.upper() for e in present})


def phase_has_any(dbf: Database, phase: str, comps):
    pobj = dbf.phases[phase]
    const = pobj.constituents
    species = {str(s).upper() for subl in const for s in subl}
    return bool(species & {c.upper() for c in comps})


def stable_phase_fractions(eq: xr.Dataset):
    """Return mapping {phase_name: fraction} at the single state point.

    Supports modern pycalphad datasets where NP and Phase are indexed by 'vertex'.
    """
    assert 'NP' in eq.data_vars, "Equilibrium dataset missing NP"
    NP = eq['NP'].squeeze()

    # Extract per-vertex phase labels if available
    phase_labels = None
    if 'Phase' in eq.data_vars:
        Pvar = eq['Phase'].squeeze()
        # decode bytes to str
        arr = np.asarray(Pvar.values)
        phase_labels = [p.decode() if isinstance(p, (bytes, bytearray)) else str(p) for p in arr]

    # At a single point, NP may have 'vertex' dim; otherwise it's scalar/empty
    amounts = {}
    if 'vertex' in NP.dims and phase_labels is not None:
        NPvals = np.asarray(NP.values)
        # Flatten along vertex
        for idx, name in enumerate(phase_labels):
            amt = float(NPvals[..., idx])
            if not np.isfinite(amt) or amt <= 0:
                continue
            amounts[name] = amounts.get(name, 0.0) + amt
    else:
        # Fallback: treat as total moles in an unknown single phase
        val = float(np.asarray(NP.values).squeeze())
        amounts['TOTAL'] = max(val, 0.0)

    tot = sum(amounts.values())
    if tot <= 0:
        return {}
    return {k: v / tot for k, v in amounts.items()}


def main():
    # Absolute path to TDB directory (from CALPHAD.md guidance)
    tdb_dir = Path('/Users/delip/play/miniclaimspy/TnETDBDB')

    # System definition
    elements = ['AL', 'CU', 'LI']
    comps = elements + ['VA']
    T = 298.0  # K
    P = 101325  # Pa
    # N components => specify N-1 mole fractions; X(AL) implied
    conds = {v.T: T, v.P: P, v.N: 1.0, v.X('CU'): 0.80, v.X('LI'): 0.10}

    # Discover and vet TDBs
    candidates = []
    considered = []
    for tdb in find_local_tdbs(tdb_dir):
        considered.append(str(tdb))
        try:
            dbf = Database(str(tdb))
        except Exception as e:
            # Skip unreadable databases
            continue
        if db_covers_elements(dbf, elements):
            candidates.append((tdb, dbf))

    if not candidates:
        msg = {
            'error': 'No suitable TDB found in TnETDBDB covering elements AL, CU, LI.',
            'tdb_dir': str(tdb_dir),
            'considered': considered,
        }
        print(json.dumps(msg, indent=2))
        sys.exit(2)

    # Choose first suitable DB by simple heuristic (only one expected: COST507R)
    chosen_path, db = candidates[0]

    # Build an inclusive phase list with any of the target elements
    phases = sorted(p for p in db.phases if phase_has_any(db, p, elements))
    if not phases:
        print(json.dumps({
            'error': 'No phases selected; did you build the phases list?',
            'chosen_db': str(chosen_path)
        }, indent=2))
        sys.exit(3)

    # Compute equilibrium at a single state point
    try:
        eq = equilibrium(db, comps, phases, conds, calc_opts={'pdens': 1500})
    except Exception as e:
        print(json.dumps({'error': f'Equilibrium calculation failed: {e}',
                          'chosen_db': str(chosen_path)}, indent=2))
        sys.exit(4)

    # Phase fractions
    # Debug: write minimal eq summary to assist parsing issues
    try:
        dbg = {
            'eq_dims': dict(eq.dims),
            'eq_coords': list(eq.coords.keys()),
            'eq_data_vars': list(eq.data_vars.keys()),
            'NP_dims': tuple(eq['NP'].dims) if 'NP' in eq.data_vars else None,
            'Phase_dims': tuple(eq['Phase'].dims) if 'Phase' in eq.data_vars else None,
        }
        Path('cu80li10al10_eq_debug.json').write_text(json.dumps(dbg, indent=2))
    except Exception:
        pass

    frac_map = stable_phase_fractions(eq)
    # Bundle and sort by descending fraction
    pairs = sorted(frac_map.items(), key=lambda t: t[1], reverse=True)
    # Filter NaNs and tiny negatives due to numerical noise
    pairs = [(n, (0.0 if not np.isfinite(a) or a < 0 else a)) for n, a in pairs]
    pairs.sort(key=lambda t: t[1], reverse=True)

    # Write CSV and JSON
    csv_path = Path('cu80li10al10_eq_298K.csv')
    json_path = Path('cu80li10al10_eq_298K.json')

    with csv_path.open('w') as f:
        f.write('phase,fraction\n')
        for n, a in pairs:
            f.write(f'{n},{a:.8f}\n')

    results = {
        'conditions': {'T_K': T, 'P_Pa': P},
        'composition_atfrac': {'CU': 0.80, 'LI': 0.10, 'AL': 0.10},
        'chosen_db': str(chosen_path),
        'considered_dbs': considered,
        'phases': [{'phase': n, 'fraction': a} for n, a in pairs if a > 1e-6],
        'phase_presence_threshold': 1e-6,
    }

    with json_path.open('w') as f:
        json.dump(results, f, indent=2)

    # Human summary
    # Consider phases with >1e-3 as present for reporting
    present = [(n, a) for n, a in pairs if a > 1e-3]
    total_present = len(present)
    classification = 'single-phase' if total_present == 1 else ('two-phase' if total_present == 2 else 'multi-phase')
    print(f"Database: {chosen_path}")
    print(f"State: T = {T:.1f} K, P = {P} Pa, x(CU,LI,AL) = (0.80, 0.10, 0.10)")
    print(f"Equilibrium: {classification}; phases (>0.1%):")
    for n, a in present:
        print(f"  - {n}: {a:.4f}")
    print(f"Wrote: {csv_path} and {json_path}")


if __name__ == '__main__':
    main()
