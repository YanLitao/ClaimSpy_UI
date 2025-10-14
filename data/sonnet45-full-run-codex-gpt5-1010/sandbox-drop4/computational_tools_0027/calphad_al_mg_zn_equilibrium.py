#!/usr/bin/env python3
# Equilibrium phase assemblage for Al-88 at.%, Mg-8 at.%, Zn-4 at.% vs temperature
# One-shot CALPHAD workflow per @CALPHAD.md: load DB -> set conditions -> solve -> summarize -> write JSON/CSV

import os
import sys
import json
import numpy as np
import pandas as pd

# Fail fast if pycalphad is unavailable
try:
    from pycalphad import Database, equilibrium, variables as v
except Exception as e:
    sys.stderr.write(f"ERROR: pycalphad import failed: {e}\n")
    sys.exit(2)

# Configuration
T_LOW = 298.15
T_HIGH = 975.0  # > Al liquidus (~933 K) to ensure crossing the liquidus
N_T = 151
PRESSURE = 101325.0
COMPOSITION = { 'AL': 0.88, 'MG': 0.08, 'ZN': 0.04 }
ELEMENTS = ['AL','MG','ZN','VA']  # include 'VA' to satisfy phase models
# Known local TDB directory discovered once (see k2code/CALPHAD.md)
TDB_DIR = '/Users/delip/play/k2code/TnETDBDB'

# Candidate phases we care about for Al–Mg–Zn (curated to avoid spurious phases)
CURATED_PHASES = [
    'LIQUID', 'FCC_A1', 'HCP_A3',  # solutions
    # Al–Mg intermetallics
    'ALMG_GAMMA', 'ALMG_BETA', 'ALMG_EPS',
    # Mg–Zn intermetallics (includes MgZn2 and related)
    'MGZN', 'MG2ZN3', 'MG7ZN3', 'MG2ZN11', 'MGZN2',
    # Laves prototypes that can host Mg/Zn
    'LAVES_C14', 'LAVES_C15', 'LAVES_C36',
]

# Outputs
RESULTS_JSON = 'results.json'
RESULTS_CSV = 'phase_fractions_vs_T.csv'
PROVENANCE_JSON = 'evidence.tdb.json'


def discover_tdbs(tdb_dir):
    paths = []
    if not os.path.isdir(tdb_dir):
        return paths
    for root, _, files in os.walk(tdb_dir):
        for fn in files:
            if fn.lower().endswith('.tdb'):
                paths.append(os.path.join(root, fn))
    return sorted(paths)


def db_has_all_elements(dbf, elems):
    # pycalphad Database.elements may be strings or Species; compare by upper-cased names
    names = set()
    for sp in dbf.elements:
        try:
            names.add(sp.name.upper())
        except AttributeError:
            names.add(str(sp).upper())
    return all(el in names for el in elems if el != 'VA')  # don't require VA explicitly


def select_phases(dbf):
    # Keep only curated phases that exist in the DB
    phases = [p for p in CURATED_PHASES if p in dbf.phases]
    # Ensure at least the essentials
    essentials = [p for p in ['LIQUID','FCC_A1'] if p in dbf.phases]
    phases = sorted(set(phases + essentials))
    if not phases:
        raise RuntimeError('No curated phases found in database for Al–Mg–Zn.')
    return phases


def run_equilibrium(tdb_path):
    dbf = Database(tdb_path)
    if not db_has_all_elements(dbf, ELEMENTS):
        raise ValueError('Database does not contain required elements (AL, MG, ZN).')
    phases = select_phases(dbf)

    # Conditions (N−1 composition constraints for N=3 components)
    conds = {
        v.T: np.linspace(T_LOW, T_HIGH, N_T),
        v.P: PRESSURE,
        v.N: 1.0,
        v.X('MG'): COMPOSITION['MG'],
        v.X('ZN'): COMPOSITION['ZN'],
    }

    # Compute
    eq = equilibrium(dbf, ELEMENTS, phases, conds, output='GM', verbose=False)

    # Phase fractions: aggregate NP over vertices grouped by Phase labels
    pf = eq['NP']
    phlab = eq['Phase']
    # Reduce singleton condition dims to get (T, vertex)
    for d in list(pf.dims):
        if d not in ('T', 'vertex'):
            pf = pf.isel({d: 0})
    for d in list(phlab.dims):
        if d not in ('T', 'vertex'):
            phlab = phlab.isel({d: 0})
    assert set(pf.dims) == {'T','vertex'} and set(phlab.dims) == {'T','vertex'}

    ts = pf['T'].values
    # Determine the set of phases that appear anywhere
    all_labels = set([str(x) for x in phlab.values.ravel().tolist() if str(x)])
    phase_names = sorted(all_labels)

    # Build per-T CSV by summing per phase label
    data_rows = []
    for iT, T in enumerate(ts):
        row = {'T_K': float(T)}
        labels = phlab.isel(T=iT).values.tolist()
        masses_raw = pf.isel(T=iT).values.tolist()
        # Treat NaNs as zero (unused vertices)
        masses = [float(m) if (m is not None and np.isfinite(m)) else 0.0 for m in masses_raw]
        tot = float(sum(masses))
        for ph in phase_names:
            s = sum(m for l, m in zip(labels, masses) if str(l) == ph)
            frac = (s / tot) if tot > 0 else 0.0
            if abs(frac) < 1e-9:
                frac = 0.0
            row[ph] = max(0.0, float(frac))
        data_rows.append(row)
    df = pd.DataFrame(data_rows)
    df.to_csv(RESULTS_CSV, index=False)

    # 298 K summary
    # Find nearest T index
    iT = int(np.argmin(np.abs(ts - T_LOW)))
    labels = phlab.isel(T=iT).values.tolist()
    masses_raw = pf.isel(T=iT).values.tolist()
    masses = [float(m) if (m is not None and np.isfinite(m)) else 0.0 for m in masses_raw]
    tot = float(sum(masses))
    summary_298 = {ph: 0.0 for ph in phase_names}
    for l, m in zip(labels, masses):
        l = str(l)
        if not l:
            continue
        summary_298[l] = summary_298.get(l, 0.0) + float(m) / tot if tot > 0 else 0.0

    # Identify primary and target intermetallics
    primary_phase = max(summary_298.items(), key=lambda kv: kv[1]) if summary_298 else (None, 0.0)

    def present(name):
        # Accept generic Laves prototypes as MgZn2 if the phase name matches or is a Laves_* phase present with nonzero fraction.
        if name.upper() == 'MGZN2':
            # Direct stoichiometric MgZn2
            return summary_298.get('MGZN2', 0.0) > 1e-3 or \
                   any(k.startswith('LAVES_') and summary_298.get(k, 0.0) > 1e-3 for k in summary_298)
        if name.upper() == 'FCC_A1':
            return summary_298.get('FCC_A1', 0.0) > 1e-3
        if name.upper() in {'AL12MG17','AL3MG2','ALMG_GAMMA'}:
            # gamma phase aliases in this DB: 'ALMG_GAMMA'
            return summary_298.get('ALMG_GAMMA', 0.0) > 1e-3
        return summary_298.get(name.upper(), 0.0) > 1e-3

    result = {
        'tdb_used': tdb_path,
        'elements': ELEMENTS,
        'phases_considered': phases,
        'temperature_range_K': [T_LOW, T_HIGH],
        'composition_x': COMPOSITION,
        'pressure_Pa': PRESSURE,
        'room_temperature_summary': {
            'T_K': T_LOW,
            'phase_fractions_mole': summary_298,
            'primary_phase': {'name': primary_phase[0], 'fraction_mole': primary_phase[1]},
            'checks': {
                'fcc_al_solid_solution_present': present('FCC_A1'),
                'mgzn2_laves_present': present('MGZN2'),
                'gamma_phase_present': present('ALMG_GAMMA'),
            }
        }
    }
    return result


def main():
    os.makedirs('.', exist_ok=True)

    # Discover TDBs and pick the Al–Mg–Zn capable one
    tdbs = discover_tdbs(TDB_DIR)
    provenance = {'tdb_dir': TDB_DIR, 'candidates': tdbs}

    selected = None
    for tdb in tdbs:
        try:
            dbf = Database(tdb)
            if db_has_all_elements(dbf, ELEMENTS):
                selected = tdb
                break
        except Exception:
            continue

    if selected is None:
        provenance['chosen'] = None
        with open(PROVENANCE_JSON, 'w') as f:
            json.dump(provenance, f, indent=2)
        result = {
            'status': 'error',
            'message': 'No suitable local TDB containing AL, MG, ZN found in TnETDBDB.',
            'provenance': provenance
        }
        with open(RESULTS_JSON, 'w') as f:
            json.dump(result, f, indent=2)
        print('ERROR: No suitable local TDB containing AL, MG, ZN found. Wrote evidence.tdb.json and results.json')
        sys.exit(1)

    provenance['chosen'] = selected
    with open(PROVENANCE_JSON, 'w') as f:
        json.dump(provenance, f, indent=2)

    try:
        result = run_equilibrium(selected)
        result['status'] = 'ok'
        with open(RESULTS_JSON, 'w') as f:
            json.dump(result, f, indent=2)
        # Print concise human summary
        rt = result['room_temperature_summary']
        phases_sorted = sorted(rt['phase_fractions_mole'].items(), key=lambda kv: kv[1], reverse=True)
        print('CALPHAD summary @ 298 K (mole fractions):')
        for name, frac in phases_sorted:
            if frac > 1e-6:
                print(f"  - {name}: {frac:.4f}")
        print(f"Primary phase: {rt['primary_phase']['name']} ({rt['primary_phase']['fraction_mole']:.4f})")
        print(f"Checks: FCC_A1={rt['checks']['fcc_al_solid_solution_present']}, MgZn2/Laves={rt['checks']['mgzn2_laves_present']}, gamma(AL-MG)={rt['checks']['gamma_phase_present']}")
        print(f"Results written: {RESULTS_JSON}, {RESULTS_CSV}, {PROVENANCE_JSON}")
    except Exception as e:
        err = {'status': 'error', 'message': str(e), 'provenance': provenance}
        with open(RESULTS_JSON, 'w') as f:
            json.dump(err, f, indent=2)
        print(f"ERROR during equilibrium: {e}")
        sys.exit(2)


if __name__ == '__main__':
    main()
