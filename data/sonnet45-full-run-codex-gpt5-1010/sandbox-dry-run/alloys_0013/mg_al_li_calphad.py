#!/usr/bin/env python3
"""
CALPHAD equilibrium for ternary Mg–Al–Li at 33.3 wt% each (1/3 by weight).

Workflow (single pass per CALPHAD.md):
  - Load local TDB (TnETDBDB/cost507R.TDB)
  - Convert wt% → mole fractions
  - Build phase list by constituents; include 'VA'
  - Set conditions with N−1 mole-fraction constraints
  - Compute equilibrium over 25–500 °C (ensure 450 °C included)
  - Compute liquidus via coarse→refined temperature bracketing if needed
  - Summarize: phases at 450 °C, liquid presence, Al–Li intermetallics, liquidus T
  - Write results.json and results.csv; print concise human summary

Units:
  - Temperature T in K (report °C where helpful)
  - Pressure P = 101325 Pa (1 bar)
  - Phase fractions are mole fractions of phase (dimensionless)

Figures of merit (FoMs):
  - Phase fraction per phase vs T
  - Presence/absence of liquid at 450 °C
  - Identification of Al–Li intermetallic phases at 450 °C
  - Liquidus temperature (K and °C) at stated composition
"""
from __future__ import annotations

import json
import csv
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import xarray as xr

from pycalphad import Database, variables as v, equilibrium


# --- Configuration (cache paths; do not re-read CALPHAD.md) ---
TDB_PATH = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB').resolve()
OUTDIR = Path('.').resolve()

# Overall composition: 33.3 wt% each of Mg, Al, Li
W_AL = 1.0 / 3.0
W_MG = 1.0 / 3.0
W_LI = 1.0 / 3.0

# Atomic weights (g/mol)
MW = {
    'AL': 26.9815385,
    'MG': 24.305,
    'LI': 6.94,
}

# Compute mole fractions from weight fractions (choose total mass = 1 g)
def weight_to_mole_fractions(w: Dict[str, float], mw: Dict[str, float]) -> Dict[str, float]:
    moles = {k: w[k] / mw[k] for k in w}
    tot = sum(moles.values())
    return {k: moles[k] / tot for k in moles}


# Temperature grids
T_MAIN_MIN, T_MAIN_MAX = 298.15, 773.15  # 25–500 °C
# Include 450 °C (723.15 K) exactly
T_450K = 450.0 + 273.15
Ts_main = np.unique(np.concatenate([
    np.arange(T_MAIN_MIN, T_MAIN_MAX + 0.1, 10.0),  # 10 K step for speed
    np.array([T_450K])
]))

# Pressure and total moles
P = 101325.0
N_TOT = 1.0

# Components: include vacancy for substitutional solution phases
COMPONENTS = ['AL', 'MG', 'LI', 'VA']


def constituents_of_phase(dbf: Database, phase: str) -> set:
    const = dbf.phases[phase].constituents
    s = set()
    for subl in const:
        for sp in subl:
            s.add(str(sp).upper())
    return s


def phases_subset_of_components(dbf: Database, components: List[str]) -> List[str]:
    comp_set = set(components)
    selected = []
    for ph in sorted(dbf.phases.keys()):
        try:
            const = constituents_of_phase(dbf, ph)
        except Exception:
            continue
        if const.issubset(comp_set):
            selected.append(ph)
    return selected


def phase_fractions(eq: xr.Dataset, phases_hint: List[str] | None = None) -> xr.DataArray:
    """Return per-phase fractions vs T (and other conditions), robust across pycalphad versions.

    The returned DataArray has dims including 'phase' and 'T'. If 'vertex' exists, it is summed out.
    """
    if 'NP' not in eq:
        raise KeyError("'NP' not found in equilibrium dataset")
    np_da = eq['NP']

    # Modern layout: explicit 'phase' dimension
    if 'phase' in np_da.dims:
        pf = np_da
        if 'vertex' in pf.dims:
            pf = pf.sum('vertex')
        tot = pf.sum('phase')
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = (pf / tot).fillna(0.0)
        return frac

    # Legacy layout: use 'Phase' labels over 'vertex'
    if 'vertex' not in np_da.dims or 'Phase' not in eq:
        raise KeyError("Legacy layout detected but missing 'vertex' or 'Phase' variable")

    phase_var = eq['Phase']
    # Determine phase names from hint or labels
    phases = list(phases_hint) if phases_hint else []
    if not phases:
        other_dims = set(phase_var.dims) - {'vertex'}
        red = phase_var
        if other_dims:
            red = red.isel(**{d: 0 for d in other_dims})
        labels = list(np.unique(np.asarray(red.values).astype(str)))
        phases = [p for p in labels if p]
    if not phases:
        raise RuntimeError("Could not infer any phase labels from dataset")

    per_phase = []
    for ph in phases:
        mask = (phase_var.astype(str) == ph)
        masked_np = np_da.where(mask, 0.0)
        per_phase.append(masked_np.sum('vertex'))

    stacked = xr.concat(per_phase, dim='phase').assign_coords({'phase': phases})
    tot = stacked.sum('phase')
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = (stacked / tot).fillna(0.0)
    return frac


def stable_mask_by_T(frac: xr.DataArray, tol: float = 1e-3) -> xr.DataArray:
    st = (frac > tol)
    other_dims = set(st.dims) - {'phase', 'T'}
    if other_dims:
        st = st.any(dim=list(other_dims))
    return st


def find_liquidus(db: Database, phases: List[str], x_al: float, x_mg: float, t_lo: float = 298.15,
                  t_hi: float = 1200.0, coarse_step: float = 25.0, tol: float = 1e-3) -> Tuple[float | None, float | None]:
    """Return (T_liquidus_K, T_liquidus_C). Uses coarse sweep then refine to 1 K."""
    cond_base = {v.P: P, v.N: N_TOT, v.X('AL'): x_al, v.X('MG'): x_mg}

    # Coarse sweep
    Ts = np.arange(t_lo, t_hi + 0.1, coarse_step)
    conds = cond_base.copy()
    conds[v.T] = Ts
    eqc = equilibrium(db, COMPONENTS, phases, conds, output=['GM'])
    fracc = phase_fractions(eqc, phases_hint=phases).mean(dim=list(set(phase_fractions(eqc, phases_hint=phases).dims) - {'phase', 'T'}))
    phases_list = [str(p) for p in fracc['phase'].values]
    liq_name = None
    for cand in ['LIQUID', 'LIQUID1', 'LIQUID_A1']:
        if cand in phases_list:
            liq_name = cand
            break
    if liq_name is None:
        # If database lacks explicit LIQUID, no liquidus can be determined here
        return None, None
    liq_frac = np.asarray(fracc.sel(phase=liq_name).values)
    idx = np.where(liq_frac > tol)[0]
    if idx.size == 0:
        # No liquid up to t_hi
        return None, None
    first_idx = int(idx[0])
    if first_idx == 0:
        T_star = float(Ts[0])
        return T_star, T_star - 273.15
    T_low = float(Ts[first_idx - 1])
    T_high = float(Ts[first_idx])

    # Refine in 1 K steps
    Ts_ref = np.arange(T_low, T_high + 0.1, 1.0)
    conds = cond_base.copy()
    conds[v.T] = Ts_ref
    eqr = equilibrium(db, COMPONENTS, phases, conds, output=['GM'])
    fracr = phase_fractions(eqr, phases_hint=phases).mean(dim=list(set(phase_fractions(eqr, phases_hint=phases).dims) - {'phase', 'T'}))
    liq_frac_r = np.asarray(fracr.sel(phase=liq_name).values)
    idx2 = np.where(liq_frac_r > tol)[0]
    if idx2.size == 0:
        # Fallback: return high bound
        return T_high, T_high - 273.15
    T_liq = float(Ts_ref[int(idx2[0])])
    return T_liq, T_liq - 273.15


def main():
    assert TDB_PATH.exists(), f"TDB not found: {TDB_PATH}"
    db = Database(str(TDB_PATH))

    # Verify database covers required elements
    db_elems = set(str(e).upper() for e in db.elements)
    required = {'AL', 'MG', 'LI'}
    missing = required - db_elems
    if missing:
        raise RuntimeError(f"Database missing elements: {sorted(missing)}")

    # Weight→mole conversion
    x = weight_to_mole_fractions({'AL': W_AL, 'MG': W_MG, 'LI': W_LI}, MW)
    X_AL = float(x['AL'])
    X_MG = float(x['MG'])
    X_LI = float(x['LI'])  # implied for constraints

    # Phase selection: include all phases to ensure solution + liquid are considered
    phases = sorted(db.phases.keys())
    assert phases, "No phases selected; did you build the phases list?"

    # Conditions for main sweep (25–500 °C)
    conds_main = {
        v.T: Ts_main,
        v.P: P,
        v.N: N_TOT,
        v.X('AL'): X_AL,
        v.X('MG'): X_MG,  # X(LI) implied
    }

    # Equilibrium calculation (single call for full temperature vector)
    eq_main = equilibrium(db, COMPONENTS, phases, conds_main, output=['GM'])
    frac_main = phase_fractions(eq_main, phases_hint=phases)

    # Reduce any extra dims for reporting
    other_dims = set(frac_main.dims) - {'phase', 'T'}
    if other_dims:
        frac_red = frac_main.mean(dim=list(other_dims))
    else:
        frac_red = frac_main

    # Phase list and liquid name detection
    phase_names = [str(p) for p in frac_red['phase'].values]
    liq_name = None
    for cand in ['LIQUID', 'LIQUID1', 'LIQUID_A1']:
        if cand in phase_names:
            liq_name = cand
            break

    # Fractions at 450 °C
    tol = 1e-3
    fr_450 = frac_red.sel(T=T_450K)
    phases_450 = []
    for ph in phase_names:
        val = float(fr_450.sel(phase=ph).values)
        if val > tol:
            phases_450.append((ph, val))
    phases_450.sort(key=lambda x: -x[1])

    # Liquid presence at 450 °C
    has_liquid_450 = False
    liquid_fraction_450 = 0.0
    if liq_name is not None:
        liquid_fraction_450 = float(fr_450.sel(phase=liq_name).values)
        has_liquid_450 = bool(liquid_fraction_450 > tol)

    # Identify Al–Li intermetallics at 450 °C (by constituents; exclude common disordered solutions and liquid)
    exclude_solution = {'FCC_A1', 'HCP_A3', 'BCC_A2', 'GAS', 'LIQUID', 'LIQUID1', 'LIQUID_A1'}
    al_li_intermetallics_450: List[Tuple[str, float]] = []
    for ph, val in phases_450:
        if ph.upper() in exclude_solution:
            continue
        try:
            const = constituents_of_phase(db, ph)
        except Exception:
            continue
        if {'AL', 'LI'}.issubset(const):
            al_li_intermetallics_450.append((ph, val))

    # Liquidus determination (may be above 500 °C)
    T_liq_K, T_liq_C = find_liquidus(db, phases, X_AL, X_MG, t_lo=298.15, t_hi=1200.0, coarse_step=25.0, tol=tol)

    # Write outputs
    json_path = OUTDIR / 'results.json'
    csv_path = OUTDIR / 'results.csv'

    # JSON
    summary = {
        'tdb_path': str(TDB_PATH),
        'components': COMPONENTS,
        'composition': {
            'wt_fraction': {'AL': W_AL, 'MG': W_MG, 'LI': W_LI},
            'mole_fraction': {'AL': X_AL, 'MG': X_MG, 'LI': X_LI},
        },
        'T_main_range_K': [float(Ts_main.min()), float(Ts_main.max())],
        'T_probe_450C_K': T_450K,
        'phases_considered': phases,
        'phases_at_450C': [{'phase': n, 'fraction': val} for n, val in phases_450],
        'has_liquid_at_450C': has_liquid_450,
        'liquid_fraction_at_450C': liquid_fraction_450,
        'al_li_intermetallics_at_450C': [{'phase': n, 'fraction': val} for n, val in al_li_intermetallics_450],
        'liquidus_T_K': T_liq_K,
        'liquidus_T_C': T_liq_C,
    }
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)

    # CSV of phase fractions vs T (reduced dims)
    phases_csv = [str(p) for p in frac_red['phase'].values]
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['T_K'] + phases_csv)
        for T in frac_red['T'].values:
            row = [float(T)] + [float(frac_red.sel(T=T, phase=ph).values) for ph in phases_csv]
            w.writerow(row)

    # Human-readable summary
    print('Database:', str(TDB_PATH))
    print('Components:', COMPONENTS)
    print('Phases considered (n=%d): %s%s' % (len(phases), ', '.join(phases[:20]), '...' if len(phases) > 20 else ''))
    print('Composition (wt frac): AL=%.3f, MG=%.3f, LI=%.3f' % (W_AL, W_MG, W_LI))
    print('Composition (mole frac): AL=%.5f, MG=%.5f, LI=%.5f' % (X_AL, X_MG, X_LI))
    print('Temperature range (main): %.2f–%.2f K (Δ10 K); probe at 450 °C' % (Ts_main.min(), Ts_main.max()))

    print('At 450 °C (723.15 K):')
    if phases_450:
        for n, val in phases_450:
            print('  %s: %.4f' % (n, val))
    else:
        print('  (no stable phases above tol=1e-3)')
    if liq_name is not None:
        print('  Liquid present? %s (fraction %.4f)' % ('Yes' if has_liquid_450 else 'No', liquid_fraction_450))
    else:
        print('  Liquid phase name not found in DB (cannot assess liquid at 450 °C).')

    if al_li_intermetallics_450:
        names = ', '.join(f"{n} ({val:.4f})" for n, val in al_li_intermetallics_450)
        print('  Al–Li intermetallics: %s' % names)
    else:
        print('  Al–Li intermetallics: none above tol=1e-3')

    if T_liq_K is not None:
        print('Estimated liquidus: %.2f K (%.2f °C)' % (T_liq_K, T_liq_C))
    else:
        print('Estimated liquidus: not reached within coarse sweep (≤ 1200 K).')

    print('Wrote:', str(json_path))
    print('Wrote:', str(csv_path))


if __name__ == '__main__':
    main()
