"""
Al–Be–W CALPHAD evaluation for the claim:

- Base alloy: AlBeMet AM162 (62 wt% Be, 38 wt% Al)
- Additions: +1, +2, +3 wt% W (final composition by weight)

Questions addressed:
  (1) Liquidus temperature of the base 62Be–38Al alloy
  (2) Whether 1–3 wt% W increases the liquidus
  (3) Phases present at 700 °C (973.15 K) in base alloy and with W
  (4) Whether W-containing phases exist at 700 °C that could pin grain boundaries

Method: Single-pass CALPHAD per CALPHAD.md guidance (load DB → conditions → solve → summarize → write JSON/CSV, no repeated DB scans). Degrees of freedom set with N−1 mole-fraction constraints.

Database: COST507R (local) located in TnETDBDB.

Outputs:
  - results.json: structured results for base and W-bearing compositions
  - results.csv: tabular summary

Notes:
  - If the database lacks a required component/interaction (e.g., Be or W), the script fails fast and reports a clear message. Base-only results are provided if W is missing but Al–Be exists; nothing is attempted if Be is missing.
  - Phase list is restricted to those whose constituents are a subset of {AL, BE, W, VA}, and must include a liquid-like phase containing the components.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pycalphad import Database, equilibrium, variables as v


# --- Configuration (cache paths; do not re-read CALPHAD.md) ---
TDB_PATH = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB').resolve()
T_TARGET = 973.15  # 700 °C in K
P_REF = 101325  # Pa


@dataclass
class CompCase:
    label: str
    w_be: float
    w_al: float
    w_w: float


def weight_to_mole_fractions(w: Dict[str, float], M: Dict[str, float]) -> Dict[str, float]:
    n = {el: w[el] / M[el] for el in w}
    tot = sum(n.values())
    return {el: (n[el] / tot) for el in n}


def pick_phases_for(dbf: Database, allowed: set) -> List[str]:
    def phase_ok(phase: str) -> bool:
        const = dbf.phases[phase].constituents
        elems = {c for subl in const for c in subl}
        return elems.issubset(allowed)

    phases = sorted([p for p in dbf.phases if phase_ok(p)])
    # Ensure at least one liquid-like phase is present
    liqs = [p for p in phases if 'LIQ' in p.upper()]
    if not liqs:
        # If filtered out, add any phase with 'LIQ' substring from the DB that contains allowed comps
        liqs = [p for p in dbf.phases if 'LIQ' in p.upper() and p not in phases]
        phases.extend(liqs)
    return phases


def get_liquidus_temperature(eq, phases: List[str]) -> float | None:
    # Phase moles NP; sum over vertex if present
    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')

    # Identify a liquid-like phase among the selected phases
    liq_candidates = [p for p in phases if 'LIQ' in p.upper() and p in pf['phase'].values]
    if not liq_candidates:
        return None
    liq = liq_candidates[0]

    # liquid fraction vs T
    phase_axis = pf['phase']
    liq_idx = int(np.where(phase_axis.values == liq)[0][0])
    liq_frac = pf.isel(phase=liq_idx)
    # Collapse remaining dims except T
    for d in list(liq_frac.dims):
        if d not in {'T',}:
            liq_frac = liq_frac.mean(d)

    Tvals = liq_frac['T'].values
    fvals = liq_frac.values

    # Scan from high T to low T; find first where solid fraction > eps
    eps = 1e-3
    solid = 1 - fvals
    # If at max T, solid already > eps, no true liquidus in scan window
    # Find first index where solid > eps
    idxs = np.where(solid > eps)[0]
    if len(idxs) == 0:
        return None
    i = idxs[0]
    if i == 0:
        return None

    # Linear interpolation between points i-1 (solid ~ 0) and i (solid > eps)
    T1, s1 = Tvals[i-1], solid[i-1]
    T2, s2 = Tvals[i], solid[i]
    if s2 == s1:
        return float(T2)
    frac = (eps - s1) / (s2 - s1)
    Tliq = T1 + frac * (T2 - T1)
    return float(Tliq)


def stable_phases_at_T(eq, Tsel: float, threshold: float = 1e-3) -> List[Tuple[str, float]]:
    # Interpolate or select nearest T index
    Tgrid = eq['T'].values
    idx = int(np.abs(Tgrid - Tsel).argmin())

    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')

    pf_T = pf.isel(T=idx)
    # Sum over any remaining non-phase dims
    for d in list(pf_T.dims):
        if d not in {'phase'}:
            pf_T = pf_T.sum(d)

    phases = []
    for j, ph in enumerate(pf_T['phase'].values):
        frac = float(pf_T.isel(phase=j).values)
        if frac > threshold:
            phases.append((str(ph), frac))
    # Sort by decreasing fraction
    phases.sort(key=lambda x: -x[1])
    return phases


def run_case(dbf: Database, case: CompCase, phases: List[str]) -> Dict:
    # Atomic weights (g/mol)
    M = {'AL': 26.9815385, 'BE': 9.0121831, 'W': 183.84}
    w = {'AL': case.w_al, 'BE': case.w_be, 'W': case.w_w}
    x = weight_to_mole_fractions(w, M)

    # Build pycalphad conditions with N−1 mole fraction constraints
    comps = ['AL', 'BE'] + (['W'] if case.w_w > 0 else [])
    conds = {
        v.T: np.linspace(2200, 600, 162),
        v.P: P_REF,
        v.N: 1.0,
    }
    # Enforce N−1 constraints in order except the last component
    # Use only components present in composition
    comps_for_constraints = comps[:-1]
    # Map to mole fractions
    for el in comps_for_constraints:
        conds[v.X(el)] = x[el]

    # Sanity check for constraints
    if len([k for k in conds if hasattr(k, 'X')]) != len(comps) - 1:
        raise ValueError('Composition constraints mismatch (N−1 rule).')

    # Equilibrium
    eq = equilibrium(dbf, ['AL','BE'] + (['W'] if case.w_w > 0 else []), phases, conds, verbose=False)

    # Liquidus
    Tliq = get_liquidus_temperature(eq, phases)

    # Stable phases at T_TARGET
    phases_700 = stable_phases_at_T(eq, T_TARGET)

    # Report W-containing phases at T_TARGET
    w_phases_700 = [(ph, frac) for ph, frac in phases_700 if 'W' in ph.upper()]

    return {
        'label': case.label,
        'mass_fractions': w,
        'mole_fractions': x,
        'T_liquidus_K': Tliq,
        'phases_700C': [{'phase': ph, 'fraction': f} for ph, f in phases_700],
        'w_phases_700C': [{'phase': ph, 'fraction': f} for ph, f in w_phases_700],
    }


def main():
    out_json = Path('results.json')
    out_csv = Path('results.csv')

    assert TDB_PATH.exists(), f"TDB not found at {TDB_PATH}"
    dbf = Database(str(TDB_PATH))

    have_al = 'AL' in dbf.elements
    have_be = 'BE' in dbf.elements
    have_w = 'W' in dbf.elements

    # Select phases limited to subset of {AL, BE, W, VA}
    allowed = {'AL','BE','W','VA'}
    phases = pick_phases_for(dbf, allowed)

    # Fail fast if Be is missing: cannot address the problem at all
    if not (have_al and have_be):
        summary = {
            'tdb_used': str(TDB_PATH),
            'db_elements': sorted(dbf.elements),
            'error': 'Database lacks AL and/or BE; cannot evaluate Al–Be system.'
        }
        out_json.write_text(json.dumps(summary, indent=2))
        # Also emit a minimal CSV for consistency
        pd.DataFrame([
            {
                'label': 'Base_62Be_38Al',
                'T_liquidus_K': np.nan,
                'phases_700C': summary['error'],
                'has_W_phase_700C': ''
            }
        ]).to_csv(out_csv, index=False)
        print(json.dumps(summary, indent=2))
        return

    cases: List[CompCase] = []
    # Base (no W)
    cases.append(CompCase(label='Base_62Be_38Al', w_be=0.62, w_al=0.38, w_w=0.0))
    # With W additions (final composition by weight; remaining split 62:38)
    for wW in (0.01, 0.02, 0.03):
        rem = 1.0 - wW
        cases.append(CompCase(label=f"{int(wW*100)}wtW", w_be=0.62*rem, w_al=0.38*rem, w_w=wW))

    results = []
    for case in cases:
        if case.w_w > 0 and not have_w:
            results.append({
                'label': case.label,
                'skipped': True,
                'reason': 'Database lacks W; cannot evaluate W additions',
            })
            continue
        try:
            res = run_case(dbf, case, phases)
        except Exception as e:
            res = {'label': case.label, 'error': str(e)}
        results.append(res)

    # Summaries
    # CSV: label, T_liquidus_K, phases_700C (semicolon-separated phase:frac), w_phase_flag
    rows = []
    for r in results:
        if r.get('skipped') or r.get('error'):
            rows.append({
                'label': r['label'],
                'T_liquidus_K': np.nan,
                'phases_700C': r.get('reason') or r.get('error'),
                'has_W_phase_700C': '',
            })
            continue
        ph700 = ';'.join([f"{p['phase']}:{p['fraction']:.4f}" for p in r['phases_700C']])
        has_w = 'yes' if len(r['w_phases_700C']) > 0 else 'no'
        rows.append({
            'label': r['label'],
            'T_liquidus_K': r['T_liquidus_K'],
            'phases_700C': ph700,
            'has_W_phase_700C': has_w,
        })

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)

    payload = {
        'tdb_used': str(TDB_PATH),
        'db_elements': sorted(dbf.elements),
        'phases_used': phases,
        'results': results,
    }
    out_json.write_text(json.dumps(payload, indent=2))

    # Human summary
    print("CALPHAD summary (COST507R):")
    print(f"TDB: {TDB_PATH}")
    print(f"Elements in DB: {', '.join(sorted(dbf.elements))}")
    print(f"Phases considered: {', '.join(phases[:15])}{' ...' if len(phases)>15 else ''}")
    for r in results:
        print(f"-- {r['label']}")
        if r.get('skipped'):
            print(f"   SKIPPED: {r['reason']}")
            continue
        if r.get('error'):
            print(f"   ERROR: {r['error']}")
            continue
        print(f"   Liquidus: {r['T_liquidus_K']:.1f} K" if r['T_liquidus_K'] else "   Liquidus: not found in scan range")
        ph700 = ', '.join([f"{p['phase']} ({p['fraction']:.3f})" for p in r['phases_700C']])
        print(f"   700°C phases: {ph700 if ph700 else 'none'}")
        if r['w_phases_700C']:
            print(f"   W-bearing phases at 700°C: {', '.join([p['phase'] for p in r['w_phases_700C']])}")


if __name__ == '__main__':
    main()
