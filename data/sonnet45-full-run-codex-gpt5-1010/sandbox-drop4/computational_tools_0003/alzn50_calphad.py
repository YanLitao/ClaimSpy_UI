#!/usr/bin/env python3
"""
Equilibrium phases for Al50Zn50 (50 at.% Al, 50 at.% Zn) from 1–700 K
using pycalphad and the local COST507R.TDB database.

Outputs:
- results.json: per-temperature phase fractions + stable phases
- results.csv: long-form table of phase fractions vs T

Notes:
- Uses absolute path to TDB under TnETDBDB, per CALPHAD.md.
- Fixes degrees of freedom: binary substitutional → one X constraint.
- Prints a concise human summary (phases at representative Ts).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from pycalphad import Database, Model, variables as v
from pycalphad import equilibrium


def phase_has_any(dbf: Database, phase: str, species: set) -> bool:
    """Return True if any constituent across sublattices intersects species."""
    const = dbf.phases[phase].constituents
    flat = set()
    for subl in const:
        flat.update(str(x) for x in subl)
    return bool(flat & species)


def main():
    # Absolute path to local TDB (do not download). Readme states scope includes Al, Zn.
    tdb_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    if not tdb_path.exists():
        print(f"ERROR: TDB not found at {tdb_path}", file=sys.stderr)
        sys.exit(2)

    # Load database once
    dbf = Database(str(tdb_path))

    # Components and allowed species
    components = ['AL', 'ZN', 'VA']
    species = set(components)

    # Start with a conservative whitelist for Al–Zn system, then intersect with DB
    preferred = ['LIQUID', 'FCC_A1', 'HCP_A3', 'HCP_ZN', 'BCC_A2']
    phases = [p for p in preferred if p in dbf.phases]
    # Fallback: if none found (unlikely), include any phase containing our species
    if not phases:
        phases = sorted(p for p in dbf.phases.keys() if phase_has_any(dbf, p, species))
    assert phases, "No phases selected; did you build the phases list?"

    # Conditions: Binary substitutional → specify N−1 mole fraction constraints
    # 50 at.% Al, 50 at.% Zn
    T_low = 1.0  # K; avoid T=0 for numerical stability
    T_high = 700.0
    T_points = np.linspace(T_low, T_high, 141)  # 5 K step approx
    conds = {v.T: T_points, v.P: 101325.0, v.N: 1.0, v.X('AL'): 0.50}

    # Compute equilibrium once over the T grid
    # Prepare outputs: long-form CSV (stable phases only) and structured JSON
    thresh = 1e-3  # phase-fraction threshold for "present"
    df_list = []
    stable_by_T = {}
    all_phase_names = set()

    # Build a temperature list that includes coarse grid + representative exact points
    reps = [300.0, 500.0, 650.0]
    coarse = np.linspace(T_low, T_high, 71)  # ~10 K step
    T_eval = np.unique(np.concatenate([coarse, np.array(reps)])).tolist()

    for T in T_eval:
        conds_T = {v.T: float(T), v.P: 101325.0, v.N: 1.0, v.X('AL'): 0.50}
        eqT = equilibrium(dbf, components, phases, conds_T, output='GM', calc_opts={'pdens': 200})

        amounts = np.array(eqT.NP.values).astype(float).squeeze()
        phases_now = eqT.Phase.values.squeeze()
        # Ensure 1D arrays over vertices
        amounts = np.atleast_1d(amounts)
        phases_now = np.atleast_1d(phases_now)

        # Map vertex-wise amounts to phase names, ignoring NaN or empty names
        phase_amount = {}
        for p_raw, val in zip(phases_now, amounts):
            name = p_raw.decode() if isinstance(p_raw, (bytes, bytearray)) else str(p_raw)
            if not name or not np.isfinite(val) or val <= 0:
                continue
            phase_amount[name] = phase_amount.get(name, 0.0) + float(val)
        all_phase_names.update(phase_amount.keys())

        tot = sum(phase_amount.values())
        # Keep only phases above threshold (by fraction)
        present = []
        if tot > 0:
            for ph, val in sorted(phase_amount.items(), key=lambda kv: kv[1], reverse=True):
                frac = val / tot
                if float(frac) > thresh:
                    present.append(ph)
                    df_list.append({'T': float(T), 'phase': ph, 'fraction': float(frac)})
        stable_by_T[str(float(T))] = present

    df = pd.DataFrame(df_list)

    rep_info = {}
    for Tr in reps:
        # Use exact T if evaluated; otherwise nearest
        if (df['T'] == Tr).any():
            Tsel = Tr
        else:
            # nearest T in df
            Ts = df['T'].unique()
            if len(Ts) == 0:
                Tsel = Tr
            else:
                Tsel = float(Ts[np.argmin(np.abs(Ts - Tr))])
        sub = df[df['T'] == float(Tsel)].copy()
        sub = sub[sub['fraction'] > thresh].sort_values('fraction', ascending=False)
        rep_info[str(Tsel)] = [{'phase': r['phase'], 'fraction': r['fraction']} for _, r in sub.iterrows()]

    # Write outputs
    out_json = {
        'tdb_path': str(tdb_path),
        'components': components,
        'phases_considered': phases,
        'phases_seen': sorted(all_phase_names),
        'T_range_K': [float(T_low), float(T_high)],
        'T_points': [float(t) for t in T_points],
        'threshold_present': thresh,
        'stable_phases_by_T': stable_by_T,
        'representative_T_summaries': rep_info,
    }
    Path('results.json').write_text(json.dumps(out_json, indent=2))
    df.to_csv('results.csv', index=False)

    # Human summary: endpoint assemblages and representative Ts
    def summarize_T(Tquery: float):
        Ts = df['T'].unique()
        if len(Ts) == 0:
            return float(Tquery), []
        Tsel = float(Tquery) if (df['T'] == float(Tquery)).any() else float(Ts[np.argmin(np.abs(Ts - Tquery))])
        sub = df[df['T'] == float(Tsel)].copy()
        sub = sub[sub['fraction'] > thresh].sort_values('fraction', ascending=False)
        return Tsel, [(r['phase'], r['fraction']) for _, r in sub.iterrows()]

    Tlo, lo_list = summarize_T(T_low)
    Thi, hi_list = summarize_T(T_high)

    print("Database:", tdb_path)
    print("Components:", components)
    print("Phases considered (subset of db):", ', '.join(phases))
    print(f"Evaluated temperatures: {len(T_eval)} points in {T_low:.1f}–{T_high:.1f} K")

    # Representative temperatures
    for Tr in reps:
        Tsel, plist = summarize_T(Tr)
        if len(plist) == 0:
            print(f"T={Tsel:.1f} K: no stable phases above {thresh}")
        else:
            assemb = ', '.join([f"{ph} ({frac:.3f})" for ph, frac in plist])
            print(f"T={Tsel:.1f} K: {assemb}")

    # Endpoints
    assemb_lo = ', '.join([f"{ph} ({frac:.3f})" for ph, frac in lo_list]) if lo_list else '—'
    assemb_hi = ', '.join([f"{ph} ({frac:.3f})" for ph, frac in hi_list]) if hi_list else '—'
    print(f"Endpoints: {Tlo:.1f} K → {assemb_lo}; {Thi:.1f} K → {assemb_hi}")

    # Single-phase vs multi-phase check across entire range (solid phases only)
    # A phase is solid if not LIQUID or GAS
    solid_names = {p for p in phases if p.upper() not in {'LIQUID', 'GAS'}}
    multi_Ts = []
    for T in T_eval:
        sub = df[(df['T'] == float(T)) & (df['phase'].isin(solid_names)) & (df['fraction'] > thresh)]
        if len(sub) > 1:
            multi_Ts.append(float(T))
    if multi_Ts:
        print(f"Two or more solid phases coexist at some T (examples near K): {multi_Ts[:5]} ...")
    else:
        print("Single solid phase across entire T range (no solid–solid two-phase detected above thresh).")


if __name__ == '__main__':
    main()
