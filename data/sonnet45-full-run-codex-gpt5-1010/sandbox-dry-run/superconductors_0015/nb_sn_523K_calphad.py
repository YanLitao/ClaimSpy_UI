#!/usr/bin/env python3
"""
Nb–Sn isothermal phase equilibria at 523 K (250 °C) using pycalphad.

Workflow:
- Load local COST507R TDB from TnETDBDB (absolute path fixed below).
- Select only phases whose constituents are a subset of {NB,SN,VA}.
- Compute equilibrium across the full composition range (X_SN ∈ [0,1]) at 523 K.
- Post-process per-phase fractions; identify stable phases across composition.
- Specifically check stability window of the A15/NB3SN phase at 523 K.
- Write results to results.json, results.csv, and a PNG plot; print concise summary.

Notes:
- Pressure fixed at 101325 Pa (1 atm).
- Phase-fraction presence threshold = 1e-3 (0.1%).
- Grid: 101 composition points for speed/coverage balance.
"""

import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dotenv import load_dotenv
from pycalphad import Database, equilibrium, variables as v


def main():
    # Environment and constants
    load_dotenv(dotenv_path=Path('.env'))  # MP_API_KEY not needed here, but loaded if present

    # Database provenance (absolute path; mandated by CALPHAD.md)
    tdb_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    assert tdb_path.exists(), f"TDB not found: {tdb_path}"

    # System setup
    T = 523.0  # K (250 °C)
    P = 101325  # Pa
    comps = ['NB', 'SN', 'VA']

    # Composition grid (binary ⇒ specify 1 mole fraction; balance implied)
    x_sn = np.linspace(0.0, 1.0, 101)
    conds = {v.T: T, v.P: P, v.N: 1.0, v.X('SN'): x_sn}

    # Load DB and select relevant phases (constituents subset of {NB,SN,VA})
    dbf = Database(str(tdb_path))
    # Quick capability check: does DB contain an A15/NB3Sn-like phase that includes Sn?
    a15_like_in_db = []
    for pname, pobj in dbf.phases.items():
        pname_u = str(pname).upper()
        species = {str(sp) for subl in pobj.constituents for sp in subl}
        if ('A15' in pname_u or 'NB3SN' in pname_u):
            a15_like_in_db.append((pname, species))

    a15_sn_supported = any(('SN' in species) for _, species in a15_like_in_db)
    if not a15_sn_supported:
        # Database not suitable for assessing Nb3Sn; emit error JSON and exit.
        error_obj = {
            'database': str(tdb_path),
            'temperature_K': T,
            'pressure_Pa': P,
            'components': comps,
            'error': 'Database lacks an A15/NB3Sn phase with Sn among constituents. Cannot assess Nb3Sn stability.',
            'a15_like_phases_found': [name for name, _ in a15_like_in_db],
        }
        with open('nb_sn_523K_error.json', 'w') as f:
            json.dump(error_obj, f, indent=2)
        print('FATAL: Selected TDB is unsuitable for Nb–Sn A15/Nb3Sn assessment.')
        print('Wrote nb_sn_523K_error.json with details.')
        return

    allowed = {'NB', 'SN', 'VA'}
    selected_phases = []
    for pname, pobj in dbf.phases.items():
        species = {str(sp) for subl in pobj.constituents for sp in subl}
        if species.issubset(allowed):
            selected_phases.append(pname)

    selected_phases = sorted(set(selected_phases))
    # If strict filtering yields nothing (DB naming peculiarities), allow solver to choose phases
    use_all_phases = False
    if not selected_phases:
        use_all_phases = True

    # Compute equilibrium (single pass)
    # Keep calc_opts modest for speed; increase pdens if convergence suffers.
    eq = equilibrium(dbf, comps, (None if use_all_phases else selected_phases), conds)

    # Phase fractions (sum vertices if present)
    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')
    tot = pf.sum('phase')
    phase_frac = (pf / tot)

    # Stable phases (threshold)
    PRES_THR = 1e-3
    stable_mask = phase_frac > PRES_THR

    # Extract axes and phase names
    phases = [p.decode() if isinstance(p, bytes) else p for p in list(phase_frac['phase'].values)]
    xvals = x_sn

    # Build a per-composition record of fractions
    # phase_frac dims typically: (T, X, P, phase) or similar ordering
    # We index at scalar T and P; keep X varying
    # Determine indices for scalar dims
    # Safely select first index of T and P
    sel = {}
    if 'T' in phase_frac.dims:
        sel['T'] = phase_frac.coords['T'].values[0]
    if 'P' in phase_frac.dims:
        sel['P'] = phase_frac.coords['P'].values[0]
    # The composition dimension for binaries is typically 'X' with a coordinate like X_SN
    # We will iterate over the X dimension directly

    # Prepare DataFrame with one column per phase
    data = {'X_SN': xvals, 'X_NB': 1.0 - xvals}
    for ph in phases:
        data[ph] = np.zeros_like(xvals)

    # Fill phase fractions per composition
    # Extract a 2D array over (X, phase)
    sliced = phase_frac
    if sel:
        sliced = phase_frac.sel(sel)

    # Ensure we have expected dims
    # After sel, dims may be ('X', 'phase') or include extra trivials; squeeze
    arr = np.squeeze(sliced.transpose(..., 'X', 'phase').values)
    # If arr is 1D because only one phase at all points, reshape
    if arr.ndim == 1:
        arr = arr.reshape((-1, 1))

    # Map phase index to name
    for j, ph in enumerate(phases):
        col = arr[:, j] if j < arr.shape[1] else np.zeros_like(xvals)
        data[ph] = col

    df = pd.DataFrame(data)

    # Determine stable phases list per composition
    stable_by_x = []
    for i, x in enumerate(xvals):
        present = [(ph, float(df.iloc[i][ph])) for ph in phases if df.iloc[i][ph] > PRES_THR]
        present.sort(key=lambda t: -t[1])
        stable_by_x.append({
            'X_SN': float(x),
            'X_NB': float(1.0 - x),
            'phases': [p for p, _ in present],
            'fractions': [a for _, a in present]
        })

    # Identify Nb3Sn (A15) phase names present in DB
    a15_like = [p for p in phases if ('A15' == p) or ('NB3SN' in p) or (p.upper().startswith('A15'))]
    a15_windows = []
    if a15_like:
        # Combine fractions for any A15-like aliases
        a15_frac = np.zeros_like(xvals)
        for p in a15_like:
            a15_frac += df[p].values
        # Stability window where fraction exceeds threshold
        mask = a15_frac > PRES_THR
        if mask.any():
            # Extract contiguous windows
            idx = np.where(mask)[0]
            # Group consecutive indices
            start = idx[0]
            prev = idx[0]
            for k in idx[1:]:
                if k == prev + 1:
                    prev = k
                else:
                    a15_windows.append((xvals[start], xvals[prev]))
                    start = k
                    prev = k
            a15_windows.append((xvals[start], xvals[prev]))

    # Write outputs
    out_json = 'nb_sn_523K_results.json'
    out_csv = 'nb_sn_523K_results.csv'
    out_png = 'nb_sn_523K_isotherm.png'

    # Save CSV (phase fractions vs composition)
    df.to_csv(out_csv, index=False)

    # Save JSON summary
    # Resolve the actual phase list that appeared in the calculation
    ph_present = sorted(set([p.decode() if isinstance(p, bytes) else p for p in np.unique(eq['Phase'].values)]))

    result = {
        'database': str(tdb_path),
        'temperature_K': T,
        'pressure_Pa': P,
        'components': comps,
        'selected_phases_request': (None if use_all_phases else selected_phases),
        'phases_in_equilibrium_results': ph_present,
        'phase_presence_threshold': PRES_THR,
        'stable_by_composition': stable_by_x,
        'a15_phase_aliases': a15_like,
        'a15_stability_windows_X_SN': [(float(lo), float(hi)) for lo, hi in a15_windows]
    }
    with open(out_json, 'w') as f:
        json.dump(result, f, indent=2)

    # Plot stacked phase fractions for the most important phases (top-N by max fraction)
    max_by_phase = {ph: float(df[ph].max()) for ph in phases}
    # Keep phases that ever exceed 1% and take the top 8 by max fraction
    keep = [ph for ph, m in max_by_phase.items() if m > 0.01]
    keep = sorted(keep, key=lambda ph: max_by_phase[ph], reverse=True)[:8]

    plt.figure(figsize=(8, 4.5))
    bottom = np.zeros_like(xvals)
    for ph in keep:
        vals = df[ph].values
        plt.fill_between(xvals, bottom, bottom + vals, step='mid', alpha=0.7, label=ph)
        bottom = bottom + vals
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('X_SN (mole fraction)')
    plt.ylabel('Phase fraction (mole)')
    plt.title('Nb–Sn, 523 K isothermal phase fractions')
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)

    # Human-readable summary
    # Endpoints
    left = stable_by_x[0]
    right = stable_by_x[-1]
    # A15/Nb3Sn status
    a15_stable = bool(a15_windows)
    a15_summary = 'stable' if a15_stable else 'not stable'
    a15_window_text = (
        ", windows in X_SN: " + "; ".join([f"[{lo:.3f}, {hi:.3f}]" for lo, hi in a15_windows])
        if a15_stable else ''
    )

    print("Database:", tdb_path)
    print("Temperature (K):", T, " Pressure (Pa):", P)
    if use_all_phases:
        print("Selected phases: None (let solver choose; filtered by components)")
    else:
        print("Selected phases (NB–SN–VA subset):", ", ".join(selected_phases))
    print(f"Left endpoint (X_SN=0): phases={left['phases']} fractions={left['fractions']}")
    print(f"Right endpoint (X_SN=1): phases={right['phases']} fractions={right['fractions']}")
    print(f"Nb3Sn (A15) phase is {a15_summary} at 523 K{a15_window_text}.")
    print("Wrote:")
    print(" -", out_json)
    print(" -", out_csv)
    print(" -", out_png)


if __name__ == '__main__':
    main()
