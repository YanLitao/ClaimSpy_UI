#!/usr/bin/env python3
"""
Al–Mg–Zn tau-phase (T-phase) equilibrium analysis at specified compositions.

Design per CALPHAD.md workflow (cached essentials):
- Use pycalphad; load DB from /Users/delip/play/miniclaimspy/TnETDBDB
- One script: load DB → set conds → solve → summarize → write JSON/CSV
- Fix DOF: N=1, set N−1 mole fractions (X(MG)=0.345, X(ZN) ∈ {0,0.02,0.04,0.06})
- Phase fractions from eq['NP'] with vertex-sum; stable if fraction > 1e-3

Outputs in CWD:
- results_al_mg_zn_tau.json
- results_al_mg_zn_tau.csv

Notes:
- Database phase naming for tau/T varies by assessment; we detect candidates
  by phase names containing 'TAU' or clearly matching Mg32(Al,Zn)49-like naming.
- If database lacks a tau/T phase, reported tau fraction is 0.0.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    from dotenv import load_dotenv  # noqa: F401
except Exception:
    # Not required for this task; only used for MP_API_KEY in other workflows
    def load_dotenv(*args, **kwargs):
        return False

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf: Database, phase: str, comps: List[str]) -> bool:
    const = dbf.phases[phase].constituents
    flat = {c for subl in const for c in subl}
    return bool(flat & set(comps))


def detect_tau_phase_names(dbf: Database, phases: List[str]) -> List[str]:
    """Return list of phase names that likely correspond to tau/T-phase.

    Heuristics:
    - Name contains 'TAU' (common)
    - Name equals 'T' or contains 'T_PHASE'/'T-PHASE' (less common; risky)
    - Name contains patterns for Mg32(Al,Zn)49, e.g., 'MG32', 'ALZN49', 'MG32ALZN49'
    The function errs on the conservative side to avoid false positives.
    """
    cands = []
    for ph in phases:
        name = ph.upper()
        if 'TAU' in name:
            cands.append(ph)
            continue
        if 'T_PHASE' in name or 'T-PHASE' in name:
            cands.append(ph)
            continue
        # Very specific stoichiometry cues
        if 'MG32' in name and ('AL' in name and 'ZN' in name):
            cands.append(ph)
            continue
        if name in {"T", "TAU_PHASE", "TAUPHASE"}:
            cands.append(ph)
            continue
    # Deduplicate while preserving order
    seen = set()
    out = []
    for x in cands:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out


def compute_equilibria(db_path: Path) -> Tuple[Dict, List[Dict]]:
    # Load environment (for completeness; not used directly here)
    load_dotenv(dotenv_path=Path('.env'))

    # Load database
    dbf = Database(str(db_path))

    comps = ['AL', 'MG', 'ZN']
    # Consider all phases from the database to avoid over-filtering
    phases = sorted(dbf.phases.keys())
    assert phases, "No phases in database; check TDB."

    tau_names = detect_tau_phase_names(dbf, phases)

    # Conditions
    T_all = np.array([298.0, 400.0, 500.0, 600.0])
    X_MG = 0.345
    X_ZN_vals = np.array([0.00, 0.02, 0.04, 0.06])
    # X(AL) implied by N−1 constraints
    conds = {v.T: T_all, v.P: 101325.0, v.N: 1.0, v.X('MG'): X_MG, v.X('ZN'): X_ZN_vals}

    # Compute equilibrium (single call for all T and Zn)
    eq = equilibrium(dbf, comps, phases, conds, verbose=False)
    # Debug: summarize dataset structure for robust fraction extraction
    try:
        print("EQ variables:", list(eq.data_vars))
        for name in eq.data_vars:
            da = eq[name]
            print(f"  {name}: dims={da.dims}")
        if 'Phase' in eq.data_vars:
            print("Phase variable dtype:", eq['Phase'].dtype)
    except Exception:
        pass

    # We'll aggregate per-phase moles using the 'vertex' decomposition and 'Phase' labels.
    # This is robust across pycalphad versions where a 'phase' dimension may be absent.

    # Prepare outputs
    results = {
        'database': str(db_path),
        'components': comps,
        'phases_considered': phases,
        'tau_candidates': tau_names,
        'conditions': {
            'T_K': T_all.tolist(),
            'P_Pa': 101325.0,
            'X_MG': float(X_MG),
            'X_ZN': X_ZN_vals.tolist(),
        },
        'threshold_stable': 1e-3,
        'data': []  # filled with per-condition summaries
    }

    rows_csv: List[Dict[str, object]] = []

    # Build mapping from index to composition/temperature
    # eq coords: typically have coords T and X_ZN (plus internal ones)
    T_coord = eq.coords['T'].values
    ZN_coord = eq.coords['X_ZN'].values

    # For each composition and temperature, record phases present and tau fraction
    for iz, xzn in enumerate(ZN_coord):
        xal = float(1.0 - X_MG - float(xzn))
        # Determine if tau present at 298 K; we'll always compute but mark the logic
        tau_present_at_298 = False
        per_comp_summary = []
        if iz == 0:
            # Debug a single condition: baseline at 298 K
            try:
                ph0 = eq['Phase'].isel(N=0, P=0, X_MG=0, T=0, X_ZN=0).values
                np0 = eq['NP'   ].isel(N=0, P=0, X_MG=0, T=0, X_ZN=0).values
                print("Debug baseline 298K: unique phases and NP moles")
                uniq = {}
                for a, b in zip(ph0, np0):
                    uniq[str(a)] = uniq.get(str(a), 0.0) + float(b)
                items = sorted(uniq.items(), key=lambda x: -x[1])
                for k, m in items:
                    print(f"  {k}: NP={m:.6g}")
                ph1 = eq['Phase'].isel(N=0, P=0, X_MG=0, T=3, X_ZN=0).values
                np1 = eq['NP'   ].isel(N=0, P=0, X_MG=0, T=3, X_ZN=0).values
                print("Debug baseline 600K: unique phases and NP moles")
                uniq1 = {}
                for a, b in zip(ph1, np1):
                    uniq1[str(a)] = uniq1.get(str(a), 0.0) + float(b)
                items1 = sorted(uniq1.items(), key=lambda x: -x[1])
                for k, m in items1:
                    print(f"  {k}: NP={m:.6g}")
            except Exception as e:
                print("Debug failure:", e)
        for it, Tval in enumerate(T_coord):
            # Extract phase labels and moles for each vertex at this condition
            ph_lab = eq['Phase'].isel(N=0, P=0, X_MG=0, T=it, X_ZN=iz)
            np_moles = eq['NP'   ].isel(N=0, P=0, X_MG=0, T=it, X_ZN=iz)
            # Convert to numpy arrays
            ph_arr = np.array(ph_lab.values)
            np_arr = np.array(np_moles.values, dtype=float)
            # Guard against shapes without 'vertex' (shouldn't happen, but be safe)
            if 'vertex' in ph_lab.dims:
                pass
            else:
                # Treat as single-vertex
                ph_arr = np.array([str(ph_arr)])
                np_arr = np.array([float(np_arr)])

            # Aggregate moles per phase name
            per_phase_moles: Dict[str, float] = {}
            for ph_name, m in zip(ph_arr, np_arr):
                key = str(ph_name).strip()
                val = float(m)
                if not key or not np.isfinite(val):
                    continue
                per_phase_moles[key] = per_phase_moles.get(key, 0.0) + val

            total_moles = float(sum(per_phase_moles.values()))
            # Convert to fractions; filter stable phases (>1e-3)
            phase_fracs_map: Dict[str, float] = {}
            if total_moles > 0:
                for k, m in per_phase_moles.items():
                    phase_fracs_map[k] = m / total_moles

            stable_items = [(p, f) for p, f in phase_fracs_map.items() if f > 1e-3]
            stable_items.sort(key=lambda x: -x[1])
            stable_phases = [p for p, _ in stable_items]
            stable_fracs = [f for _, f in stable_items]
            total_frac = float(sum(stable_fracs))

            # Tau/T fraction
            tau_frac = 0.0
            tau_name_hit = None
            for cand in tau_names:
                if cand in phase_fracs_map:
                    tau_frac += phase_fracs_map[cand]
                    if tau_name_hit is None:
                        tau_name_hit = cand
            if Tval == 298.0 and tau_frac > 1e-9:
                tau_present_at_298 = True

            # Record per-condition
            per_comp_summary.append({
                'T_K': float(Tval),
                'X_AL': xal,
                'X_MG': float(X_MG),
                'X_ZN': float(xzn),
                'phases_present': stable_phases,
                'phase_fractions': {p: f for p, f in zip(stable_phases, stable_fracs)},
                'tau_phase_name': tau_name_hit,
                'tau_phase_fraction': tau_frac,
                'total_phase_fraction_sum': total_frac,
            })

        # Apply the reporting rule: always include 298 K; if tau absent at 298 K, include 400/500/600 K
        filtered = []
        for item in per_comp_summary:
            if item['T_K'] == 298.0:
                filtered.append(item)
            else:
                if not tau_present_at_298:
                    filtered.append(item)

        results['data'].append({
            'composition': {'X_AL': xal, 'X_MG': float(X_MG), 'X_ZN': float(xzn)},
            'tau_present_at_298K': tau_present_at_298,
            'entries': filtered,
        })

        # CSV rows for all computed temps (not just filtered) to enable easy scanning
        for item in per_comp_summary:
            for p, f in item['phase_fractions'].items():
                rows_csv.append({
                    'T_K': item['T_K'],
                    'X_AL': item['X_AL'],
                    'X_MG': item['X_MG'],
                    'X_ZN': item['X_ZN'],
                    'phase': p,
                    'fraction': f,
                    'is_tau': (p == item['tau_phase_name']),
                    'sum_fractions_this_condition': item['total_phase_fraction_sum'],
                })

    return results, rows_csv


def main():
    # Resolve TDB path (cached from prior runs and CALPHAD.md):
    TDB_PATH = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB').resolve()
    if not TDB_PATH.exists():
        raise FileNotFoundError(f"TDB not found at {TDB_PATH}")

    results, rows_csv = compute_equilibria(TDB_PATH)

    # Write outputs
    out_json = Path('results_al_mg_zn_tau.json')
    out_csv = Path('results_al_mg_zn_tau.csv')

    with out_json.open('w') as f:
        json.dump(results, f, indent=2)

    # CSV header
    import csv
    fieldnames = ['T_K', 'X_AL', 'X_MG', 'X_ZN', 'phase', 'fraction', 'is_tau', 'sum_fractions_this_condition']
    with out_csv.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_csv)

    # Human summary: phases per condition and tau fraction
    print("Database:", results['database'])
    print("Tau-phase candidates found:", results['tau_candidates'] or '(none)')
    print("Compositions (X_AL, X_MG, X_ZN) at T=298K, and if tau absent also T=400/500/600K:")
    for comp_block in results['data']:
        comp = comp_block['composition']
        print(f"- Composition: X_AL={comp['X_AL']:.4f}, X_MG={comp['X_MG']:.4f}, X_ZN={comp['X_ZN']:.4f}")
        for ent in comp_block['entries']:
            tau_pct = 100.0 * ent['tau_phase_fraction']
            tot_pct = 100.0 * ent['total_phase_fraction_sum']
            phases = ', '.join(ent['phases_present']) if ent['phases_present'] else '(none)'
            print(f"  T={ent['T_K']:.0f} K | tau={tau_pct:.2f}% | sum={tot_pct:.2f}% | phases: {phases}")
    print(f"\nWrote {out_json} and {out_csv}")


if __name__ == '__main__':
    np.random.seed(42)
    main()
