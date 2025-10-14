#!/usr/bin/env python3
"""
CALPHAD analysis for TA15 (Ti-6.5Al-2Zr-1Mo-1V in wt%)

Computes beta-phase fraction at 200 C (473 K), 1000 C (1273 K), and at the
beta-transus; also estimates the liquidus (melting onset) temperature.

Outputs:
- results.csv: sampled temperatures, per-phase fractions, including beta and liquid
- results.json: key results (fractions at 473/1273 K, beta-transus, liquidus)

Notes:
- Uses the COST507R.TDB in the TnETDBDB directory.
- Filters phases to those that include any of the components {TI, AL, ZR, MO, V}.
- Treats Ti beta as a BCC-type solution; attempts to map common names (e.g., BCC_A2, BETA).
"""

import json
import math
import os
from pathlib import Path

import numpy as np

try:
    # pycalphad stack
    from pycalphad import Database, variables as v, equilibrium
    from pycalphad.plot.utils import phase_legend
except Exception as e:
    raise SystemExit(f"pycalphad import failed: {e}\nEnsure you're using the 'readertools' environment.")


def wt_to_molfrac(wt):
    """Convert weight percent dict to mole fraction dict.

    wt: dict of element symbol -> wt% (summing to ~100)
    Returns: dict of element -> mole fraction (summing to 1)
    """
    # Atomic weights (g/mol): IUPAC 2013/2017 conventional values
    AW = {
        'TI': 47.867,
        'AL': 26.9815385,
        'ZR': 91.224,
        'MO': 95.95,
        'V': 50.9415,
    }
    # moles in 100 g alloy
    moles = {}
    total_moles = 0.0
    for el, w in wt.items():
        aw = AW[el]
        n = float(w) / aw
        moles[el] = n
        total_moles += n
    return {el: n / total_moles for el, n in moles.items()}


def phase_has_any(dbf, phase, comps):
    const = dbf.phases[phase].constituents
    # Flatten constituents across sublattices
    s = {c for subl in const for c in subl}
    return bool(s & set(comps))


def detect_beta_phase(phases):
    """Pick the database phase name to represent Ti beta.

    Preference order: explicit 'BETA'/'BETA_TI', then any 'BCC_A2' or containing 'BCC'.
    Returns a phase name or None if not found.
    """
    cand = None
    up = [p.upper() for p in phases]
    # Explicit names first
    for name in ['BETA_TI', 'BETA-TI', 'BETA', 'TI_BETA']:
        if name in up:
            return phases[up.index(name)]
    # Common BCC solution model names
    for name in ['BCC_A2', 'BCC_B2']:
        if name in up:
            return phases[up.index(name)]
    # Any phase containing BCC
    for i, p in enumerate(up):
        if 'BCC' in p:
            cand = phases[i]
            break
    return cand


def detect_alpha_phase(phases):
    """Pick the database phase name to represent Ti alpha.

    Preference: HCP_A3, ALPHA_TI, ALPHA; else any 'HCP'.
    """
    up = [p.upper() for p in phases]
    for name in ['HCP_A3', 'ALPHA_TI', 'ALPHA-TI', 'ALPHA']:
        if name in up:
            return phases[up.index(name)]
    for i, p in enumerate(up):
        if 'HCP' in p:
            return phases[i]
    return None


def main():
    # Fixed inputs
    tdb_path = '/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB'
    assert os.path.exists(tdb_path), f"TDB not found at {tdb_path}"

    # TA15 nominal wt%: Ti-6.5Al-2Zr-1Mo-1V -> Ti is balance (89.5)
    wt = {'AL': 6.5, 'ZR': 2.0, 'MO': 1.0, 'V': 1.0, 'TI': 89.5}
    x = wt_to_molfrac(wt)

    # Components (include vacancy)
    components = ['TI', 'AL', 'ZR', 'MO', 'V', 'VA']

    # Load DB and determine relevant phases
    dbf = Database(tdb_path)
    all_phases = sorted(dbf.phases.keys())
    phases = sorted(p for p in dbf.phases if phase_has_any(dbf, p, ['TI', 'AL', 'ZR', 'MO', 'V']))
    if not phases:
        # Write a helpful JSON and exit cleanly if the DB does not support the system
        summary = {
            'composition_wt_percent': wt,
            'composition_mole_fraction': {k: float(v) for k, v in x.items()},
            'database': Path(tdb_path).name,
            'phases_considered': [],
            'error': 'Selected database does not contain relevant phases for TI-AL-ZR-MO-V system. Cannot perform equilibrium calculation.'
        }
        Path('results.json').write_text(json.dumps(summary, indent=2))
        print("No relevant phases found in database for TI-AL-ZR-MO-V. Wrote results.json with diagnostic info.")
        return

    # Identify key phase names for summary
    beta_name = detect_beta_phase(phases)
    alpha_name = detect_alpha_phase(phases)

    # Temperature range for scans (K)
    T_low = 400.0
    T_high = 2200.0
    nT = 182  # ~10 K steps
    T_grid = np.linspace(T_low, T_high, nT)

    # Build conditions: set N, P, and N-1 mole fraction constraints
    # For 5 substitutional components, fix 4 X's; TI is implied.
    conds = {
        v.T: T_grid,
        v.P: 101325.0,
        v.N: 1.0,
        v.X('AL'): x['AL'],
        v.X('ZR'): x['ZR'],
        v.X('MO'): x['MO'],
        v.X('V'): x['V'],
    }

    # Solve equilibrium
    eq = equilibrium(dbf, components, phases, conds, verbose=False, broadcast=True)

    # Phase fractions (NP is phase moles per condition)
    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')
    tot = pf.sum('phase')
    phase_frac = (pf / tot)

    # Helper to get array by phase name safely
    def frac_of(phase_name):
        if phase_name is None:
            return np.zeros_like(phase_frac.isel(phase=0))
        phases_arr = phase_frac['phase'].values.tolist()
        if phase_name not in phases_arr:
            return np.zeros_like(phase_frac.isel(phase=0))
        sel = phase_frac.sel(phase=phase_name)
        # Drop extraneous dims to (T)
        for d in list(sel.dims):
            if d not in ('T',):
                sel = sel.isel({d: 0})
        return np.asarray(sel.values)

    # Extract key phase fractions across T
    T_vals = np.asarray(phase_frac['T'].values)
    beta_frac = frac_of(beta_name)
    liquid_name = None
    for cand in ['LIQUID', 'LIQ', 'LIQUID0']:
        if cand in phase_frac['phase'].values:
            liquid_name = cand
            break
    liquid_frac = frac_of(liquid_name)

    # Alpha (for transus reference)
    alpha_frac = frac_of(alpha_name)

    # Values at specific temperatures
    def interp_at(Tq, y):
        # simple nearest neighbor to avoid extrap artifacts
        idx = int(np.abs(T_vals - Tq).argmin())
        return float(y[idx])

    beta_473K = interp_at(473.0, beta_frac)
    beta_1273K = interp_at(1273.0, beta_frac)

    # Beta-transus: first T where beta is essentially single-phase
    # Define as T at which (1 - frac_beta) < 1e-3 and liquid ~0.
    others_frac = 1.0 - beta_frac - liquid_frac
    mask_single_beta = (beta_frac > 1 - 1e-3) & (liquid_frac < 1e-3)
    T_beta_transus = None
    if np.any(mask_single_beta):
        T_beta_transus = float(T_vals[mask_single_beta][0])
    else:
        # fallback: temperature at which alpha drops below 1e-3 (if available)
        if alpha_name is not None and alpha_frac.size == T_vals.size:
            idxs = np.where(alpha_frac < 1e-3)[0]
            if idxs.size:
                T_beta_transus = float(T_vals[idxs[0]])

    # Liquidus: first T where liquid fraction appears (>1e-3)
    T_liquidus = None
    if liquid_frac.size == T_vals.size:
        idxs = np.where(liquid_frac > 1e-3)[0]
        if idxs.size:
            T_liquidus = float(T_vals[idxs[0]])

    # Write CSV of phase fractions vs T
    # Collect all phase names included in the result for compact CSV
    phase_names = [str(p) for p in phase_frac['phase'].values.tolist()]
    # Construct array (nT x nphases)
    data = {}
    data['T'] = T_vals.tolist()
    for pname in phase_names:
        arr = frac_of(pname)
        data[pname] = [float(a) for a in arr]

    csv_path = Path('results.csv')
    with csv_path.open('w') as f:
        # Header with T then phases
        hdr = ['T'] + phase_names
        f.write(','.join(hdr) + '\n')
        for i in range(len(T_vals)):
            row = [f"{T_vals[i]:.6f}"] + [f"{data[p][i]:.8f}" for p in phase_names]
            f.write(','.join(row) + '\n')

    # JSON summary
    summary = {
        'composition_wt_percent': wt,
        'composition_mole_fraction': {k: float(v) for k, v in x.items()},
        'database': Path(tdb_path).name,
        'phases_considered': phases,
        'beta_phase_name': beta_name,
        'alpha_phase_name': alpha_name,
        'beta_fraction_at_473K': beta_473K,
        'beta_fraction_at_1273K': beta_1273K,
        'beta_transus_temperature_K': T_beta_transus,
        'liquidus_temperature_K': T_liquidus,
    }
    json_path = Path('results.json')
    json_path.write_text(json.dumps(summary, indent=2))

    # Human-readable summary
    def fmt(x):
        return 'n/a' if x is None else f"{x:.1f}"

    print("TA15 CALPHAD summary (COST507R)")
    print(f"- Beta phase key: {beta_name}")
    print(f"- Alpha phase key: {alpha_name}")
    print(f"- Beta fraction at 473 K: {beta_473K:.4f}")
    print(f"- Beta fraction at 1273 K: {beta_1273K:.4f}")
    print(f"- Beta-transus (K): {fmt(T_beta_transus)}")
    print(f"- Liquidus (K): {fmt(T_liquidus)}")
    print(f"Wrote: {csv_path} and {json_path}")


if __name__ == '__main__':
    main()
