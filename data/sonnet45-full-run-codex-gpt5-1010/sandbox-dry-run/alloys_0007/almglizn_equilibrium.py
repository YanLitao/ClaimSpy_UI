#!/usr/bin/env python3
"""
Equilibrium phase analysis for Al55Mg35Li5Zn5 (at.%) from 300–800 K at 1 bar.
- Database: COST507 (local) from TnETDBDB
- Method: pycalphad equilibrium with fixed overall composition (N−1 X constraints)
- Outputs: results.json, results.csv, and printed summary
"""
import json
import csv
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from pycalphad import Database, variables as v, equilibrium

# --- Config (cache paths; do not re-read CALPHAD.md) ---
TDB_PATH = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB').resolve()
OUTDIR = Path('/Users/delip/play/miniclaimspy/sandbox/alloys_0007').resolve()

# Composition (atomic/mole fractions)
X_AL = 0.55
X_MG = 0.35
X_LI = 0.05
# X_ZN implied by N-1 rule
X_ZN = 1.0 - (X_AL + X_MG + X_LI)
if abs(X_ZN - 0.05) > 1e-9:
    print(f"Warning: implied X_ZN={X_ZN:.6f} differs from expected 0.05")

# Temperature grid
T_MIN, T_MAX, T_STEP = 300.0, 800.0, 10.0
Ts = np.arange(T_MIN, T_MAX + 0.1, T_STEP)

# Pressure, moles
P = 101325.0
N_TOT = 1.0

# Components: include vacancy for substitutional phases
COMPONENTS = ['AL', 'MG', 'LI', 'ZN', 'VA']


def constituents_of_phase(dbf, phase):
    const = dbf.phases[phase].constituents
    s = set()
    for subl in const:
        for sp in subl:
            s.add(str(sp).upper())
    return s


def phases_subset_of_components(dbf, components):
    comp_set = set(components)
    selected = []
    for ph in sorted(dbf.phases.keys()):
        try:
            const = constituents_of_phase(dbf, ph)
        except Exception:
            continue
        # Include if all constituents are within components
        if const.issubset(comp_set):
            selected.append(ph)
    return selected


def compute_equilibrium():
    assert TDB_PATH.exists(), f"TDB not found: {TDB_PATH}"
    db = Database(str(TDB_PATH))

    phases = phases_subset_of_components(db, COMPONENTS)
    # If overly strict filtering yields very few phases, fall back to using all phases
    if len(phases) < 3:
        phases = sorted(db.phases.keys())
    if not phases:
        raise RuntimeError("No phases selectable with the given components; check database coverage.")

    conds = {
        v.T: Ts,
        v.P: P,
        v.N: N_TOT,
        v.X('AL'): X_AL,
        v.X('MG'): X_MG,
        v.X('LI'): X_LI,  # X(ZN) implied
    }

    # Use default solver to avoid version-mismatch issues
    eq = equilibrium(db, COMPONENTS, phases, conds, output=['GM'])
    # Debug: ensure expected structure
    # Print lightweight diagnostics (not too verbose)
    try:
        dims = {k: list(v.dims) for k, v in eq.data_vars.items()}
        print(f"EQ keys: {list(eq.data_vars.keys())}")
        print(f"Dims: {dims}")
        print(f"Selected phases (n={len(phases)}): {', '.join(phases[:20])}{'...' if len(phases)>20 else ''}")
    except Exception:
        pass
    return db, phases, eq


def phase_fractions(eq: xr.Dataset, phases_hint=None) -> xr.DataArray:
    """Return per-phase fractions vs T (and other conditions).

    Supports older pycalphad layouts where 'NP' lacks a 'phase' dim
    and the phase label is stored in the 'Phase' variable over the
    'vertex' dimension.
    """
    if 'NP' not in eq:
        raise KeyError("'NP' not found in equilibrium dataset")

    np_da = eq['NP']

    # Case 1: modern layout with explicit 'phase' dim
    if 'phase' in np_da.dims:
        pf = np_da
        if 'vertex' in pf.dims:
            pf = pf.sum('vertex')
        tot = pf.sum('phase')
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = (pf / tot).fillna(0.0)
        return frac

    # Case 2: legacy layout: group NP by Phase across 'vertex'
    if 'vertex' not in np_da.dims or 'Phase' not in eq:
        raise KeyError("Legacy layout detected but missing 'vertex' or 'Phase' variable")

    phase_var = eq['Phase']  # labels per vertex
    # Reduce any dims other than vertex and conditions so broadcasting works
    # The dims for both are the same: ('N','P','T','X_*', 'vertex')
    # Build per-phase arrays by masking NP where Phase==phase_name
    # Determine list of phase names from data (not relying solely on hint)
    try:
        # Unique labels over 'vertex' (collapse other dims first)
        red = phase_var
        other_dims = set(phase_var.dims) - {'vertex'}
        if other_dims:
            # take first along other dims to collect labels generically
            red = red.isel(**{d: 0 for d in other_dims})
        labels = list(np.unique(np.asarray(red.values).astype(str)))
    except Exception:
        labels = []

    # If hint provided, prefer its order; otherwise use labels discovered
    phases = list(phases_hint) if phases_hint else labels
    # Filter out empty/None labels
    phases = [p for p in phases if p and isinstance(p, str)]
    if not phases:
        raise RuntimeError("No phase labels could be inferred from dataset")

    # Sum moles per phase across 'vertex'
    # Build a list of DataArrays (one per phase)
    per_phase = []
    for ph in phases:
        mask = (phase_var.astype(str) == ph)
        masked_np = np_da.where(mask, 0.0)
        per_phase.append(masked_np.sum('vertex'))

    # Stack into a new DataArray with a 'phase' dimension
    stacked = xr.concat(per_phase, dim='phase')
    stacked = stacked.assign_coords({'phase': phases})

    tot = stacked.sum('phase')
    with np.errstate(divide='ignore', invalid='ignore'):
        frac = (stacked / tot).fillna(0.0)
    return frac


def stable_phases_by_T(frac: xr.DataArray, tol: float = 1e-3):
    # Boolean mask for stability
    stable = (frac > tol)
    # Reduce over all non-phase, non-T dimensions
    other_dims = set(stable.dims) - {'phase', 'T'}
    if other_dims:
        stable = stable.any(dim=list(other_dims))
    return stable


def summarize(eq: xr.Dataset, phases):
    frac = phase_fractions(eq, phases_hint=phases)
    stable = stable_phases_by_T(frac)

    Ts = eq['T'].values
    phases = list(frac['phase'].values)

    # Determine single-phase regions (exactly one phase above tol)
    tol = 1e-3
    single_mask = (frac > tol).sum(dim='phase') == 1
    # Collapse any remaining dims except T
    other_dims = set(single_mask.dims) - {'T'}
    if other_dims:
        single_mask = single_mask.any(dim=list(other_dims))

    # Identify which phase is the single phase at each T
    # For selection, reduce other dims first
    frac_reduced = frac
    other_dims = set(frac_reduced.dims) - {'phase', 'T'}
    if other_dims:
        frac_reduced = frac_reduced.mean(dim=list(other_dims))

    single_phase_name = []
    for i, T in enumerate(Ts):
        if bool(single_mask.sel(T=T).values):
            # argmax among phases
            j = int(np.nanargmax(frac_reduced.sel(T=T).values))
            single_phase_name.append(phases[j])
        else:
            single_phase_name.append(None)

    # Build contiguous intervals of single-phase behavior
    intervals = []
    start = None
    current_phase = None
    for idx, (T, sp) in enumerate(zip(Ts, single_phase_name)):
        if sp is not None:
            if start is None:
                start = T
                current_phase = sp
            elif sp != current_phase:
                # phase changed; close previous
                intervals.append({'T_min': float(start), 'T_max': float(Ts[idx-1]), 'phase': current_phase})
                start = T
                current_phase = sp
        else:
            if start is not None:
                intervals.append({'T_min': float(start), 'T_max': float(Ts[idx-1]), 'phase': current_phase})
                start = None
                current_phase = None
    if start is not None:
        intervals.append({'T_min': float(start), 'T_max': float(Ts[-1]), 'phase': current_phase})

    # Endpoint summaries
    def endpoint(Tval):
        mask = stable.sel(T=Tval)
        stables = [str(phases[i]) for i, ok in enumerate(mask.values) if ok]
        fr = frac.sel(T=Tval)
        fr_list = []
        for i, ph in enumerate(phases):
            val = float(fr.sel(phase=ph).values)
            if val > tol:
                fr_list.append((str(ph), val))
        fr_list.sort(key=lambda x: -x[1])
        return stables, fr_list

    st_300, fr_300 = endpoint(Ts[0])
    st_800, fr_800 = endpoint(Ts[-1])

    return {
        'Ts': Ts.tolist(),
        'single_phase_intervals': intervals,
        'endpoint_300K': {'stable_phases': st_300, 'fractions': [(n, float(v)) for n, v in fr_300]},
        'endpoint_800K': {'stable_phases': st_800, 'fractions': [(n, float(v)) for n, v in fr_800]},
    }, frac


def write_outputs(summary: dict, frac: xr.DataArray):
    # JSON summary
    json_path = OUTDIR / 'results.json'
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)

    # CSV of phase fractions vs T
    csv_path = OUTDIR / 'results.csv'
    phases = [str(p) for p in frac['phase'].values]
    Ts = frac['T'].values
    # reduce other dims
    other_dims = set(frac.dims) - {'phase', 'T'}
    fr = frac
    if other_dims:
        fr = fr.mean(dim=list(other_dims))

    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['T_K'] + phases)
        for T in Ts:
            row = [float(T)] + [float(fr.sel(T=T, phase=ph).values) for ph in phases]
            w.writerow(row)

    return str(json_path), str(csv_path)


def main():
    db, phases, eq = compute_equilibrium()
    summary, frac = summarize(eq, phases)

    json_path, csv_path = write_outputs(summary, frac)

    # Human summary
    print('Database:', str(TDB_PATH))
    print('Components:', COMPONENTS)
    print('Phases considered:', ', '.join(phases))
    print('Composition (mole fraction): X(AL)=%.3f, X(MG)=%.3f, X(LI)=%.3f, X(ZN)=%.3f' % (X_AL, X_MG, X_LI, X_ZN))
    print('Temperature range: %.1f–%.1f K (step %.1f K)' % (T_MIN, T_MAX, T_STEP))

    intervals = summary['single_phase_intervals']
    if intervals:
        print('Single-phase regions (tol=1e-3):')
        for itv in intervals:
            print('  %s: %.1f–%.1f K' % (itv['phase'], itv['T_min'], itv['T_max']))
    else:
        print('No single-phase region identified (tol=1e-3).')

    ep300 = summary['endpoint_300K']
    ep800 = summary['endpoint_800K']

    def fmt_ep(tag, ep):
        print(f"Stable phases at {tag}: {', '.join(ep['stable_phases']) if ep['stable_phases'] else '(none)'}")
        for name, val in ep['fractions']:
            print(f"  {name}: {val:.4f}")

    fmt_ep('300 K', ep300)
    fmt_ep('800 K', ep800)

    print('Wrote:', json_path)
    print('Wrote:', csv_path)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print('ERROR:', e)
        sys.exit(1)
