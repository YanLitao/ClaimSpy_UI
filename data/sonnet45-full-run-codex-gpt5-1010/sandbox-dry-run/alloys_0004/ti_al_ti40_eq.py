#!/usr/bin/env python3
"""
Ti–Al CALPHAD at Ti-40 at% (Al-60 at%):
 - Database: COST507R (path /Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB)
 - T range: 300–1500 K
 - P: 1 atm
 - Composition constraints: X(TI)=0.40 (X(AL)=0.60 implied)

One-shot script per CALPHAD.md guidance: load DB → set conditions → solve → summarize → write JSON/CSV.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from pycalphad import Database, equilibrium, variables as v


def phase_constituents(dbf, phase):
    const = dbf.phases[phase].constituents
    return [tuple(subl) for subl in const]


def phase_has_any(dbf, phase, comps):
    const = dbf.phases[phase].constituents
    return any({c for subl in const for c in subl} & set(comps))


def contiguous_true_ranges(mask: np.ndarray, axis_vals: np.ndarray):
    """Given boolean mask vs axis_vals (1D), return list of (start, end) inclusive ranges
    where mask is True. axis_vals must be monotonically increasing.
    """
    ranges = []
    if mask.size == 0:
        return ranges
    in_run = False
    start = None
    for i, val in enumerate(axis_vals):
        if mask[i] and not in_run:
            in_run = True
            start = val
        elif (not mask[i]) and in_run:
            in_run = False
            end = axis_vals[i - 1]
            ranges.append((float(start), float(end)))
    if in_run:
        ranges.append((float(start), float(axis_vals[-1])))
    return ranges


def main():
    # --- Inputs ---
    tdb_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    assert tdb_path.exists(), f"TDB not found: {tdb_path}"

    components = ['TI', 'AL', 'VA']
    # N-1 mole fraction constraints for an N-component (substitutional) system
    x_ti = 0.40  # Ti-40 at%, Al-60 at% implied
    p_pa = 101325  # 1 atm
    t_lo, t_hi, t_n = 300.0, 1500.0, 61  # inclusive grid
    t_grid = np.linspace(t_lo, t_hi, t_n)

    # --- Load DB and select phases ---
    dbf = Database(str(tdb_path))
    # Use all phases defined in the database; pycalphad will restrict
    # constituents to our component set automatically.
    phases = sorted(dbf.phases.keys())

    # --- Equilibrium ---
    conds = {v.T: t_grid, v.P: p_pa, v.N: 1.0, v.X('TI'): x_ti}
    eq = equilibrium(dbf, components, phases, conds, verbose=False)
    # Brief dataset structure (for robustness across pycalphad versions)
    try:
        print('Eq dims:', dict(eq.dims))
        print('Eq data vars:', list(eq.data_vars))
        print('Eq coords:', list(eq.coords))
        try:
            print('NP dims:', tuple(eq['NP'].dims))
        except Exception:
            pass
        try:
            print('Phase dims:', tuple(eq['Phase'].dims))
        except Exception:
            pass
    except Exception:
        pass

    # --- Phase fractions (robust to pycalphad layout) ---
    pf = eq['NP']              # dims: ('N','P','T','X_TI','vertex') here
    ph = eq['Phase']           # same dims, labels of phases per vertex
    # Unique phase labels observed
    labels = sorted({str(x) for x in np.unique(ph.values)})
    # Total moles across all vertices
    tot = pf.sum('vertex')     # dims: ('N','P','T','X_TI')
    # Build per-phase moles by masking NP by Phase label then summing over vertex
    moles_by_phase = {}
    for lbl in labels:
        moles_by_phase[lbl] = pf.where(ph == lbl).sum('vertex')  # dims: ('N','P','T','X_TI')
    # Prepare temperature array and phase fractions table
    Tvals = np.asarray(eq['T'])
    phase_names = labels
    frac_table = np.zeros((Tvals.size, len(phase_names)))
    for j, lbl in enumerate(phase_names):
        # squeeze to (T,) assuming singleton N,P,X_TI
        num = np.squeeze(moles_by_phase[lbl].values)
        den = np.squeeze(tot.values)
        fr = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=(den != 0))
        # Ensure shape aligns to T
        fr = fr.reshape(-1)  # T axis
        frac_table[:, j] = fr
    # Convert to xarray-like accessors via lightweight shims
    # We'll keep using numpy arrays for speed; map to previous variable names
    # phase_frac[T_index, phase_index]
    stable_mask = frac_table > 1e-3
    present_any = []
    for pi, pname in enumerate(phase_names):
        mask = stable_mask[:, pi]
        if mask.any():
            present_any.append(pname)

    # Stability ranges per phase
    stability_ranges = {}
    for pi, pname in enumerate(phase_names):
        mask = stable_mask[:, pi]
        if mask.any():
            ranges = contiguous_true_ranges(mask, Tvals)
            stability_ranges[pname] = ranges

    # Per-T summary of stable phases and phase fractions (only stable ones)
    summary_by_T = []
    for ti, T in enumerate(Tvals):
        row = {
            'T_K': float(T),
            'stable_phases': [],
            'fractions': {},
        }
        for pi, pname in enumerate(phase_names):
            frac = float(frac_table[ti, pi])
            if frac > 1e-3:
                row['stable_phases'].append(pname)
                row['fractions'][pname] = frac
        summary_by_T.append(row)

    # Heuristic detection of Ti2Al3-like phases by name
    import re
    pattern = re.compile(r'(TI\s*2\s*AL\s*3|AL\s*3\s*TI\s*2|TI\s*3\s*AL\s*2|AL\s*2\s*TI\s*3|TIAL3|TI3AL)', re.IGNORECASE)
    ti2al3_like = [p for p in phases if pattern.search(p)]

    # Endpoint summaries
    def endpoint_info(idx):
        stabs = [phase_names[pi] for pi in range(len(phase_names)) if float(frac_table[idx, pi]) > 1e-3]
        fracs = {phase_names[pi]: float(frac_table[idx, pi]) for pi in range(len(phase_names)) if float(frac_table[idx, pi]) > 1e-3}
        return stabs, fracs

    stable_300K, fr_300K = endpoint_info(0)
    stable_1500K, fr_1500K = endpoint_info(len(Tvals) - 1)

    # CSV of phase fractions vs T for phases that appear anywhere
    all_present = sorted(set(present_any))
    csv_lines = []
    header = ['T_K'] + all_present
    csv_lines.append(','.join(header))
    for ti, T in enumerate(Tvals):
        vals = [f"{T:.6g}"]
        for pname in all_present:
            pi = phase_names.index(pname)
            vals.append(f"{float(frac_table[ti, pi]):.6g}")
        csv_lines.append(','.join(vals))

    out_csv = Path('ti_al_ti40_results.csv')
    out_json = Path('ti_al_ti40_results.json')
    out_csv.write_text('\n'.join(csv_lines))

    results = {
        'tdb_used': str(tdb_path),
        'components': components,
        'composition': {'X(TI)': x_ti, 'X(AL)_implied': 1 - x_ti},
        'pressure_Pa': p_pa,
        'T_grid_K': {'T_min': float(t_lo), 'T_max': float(t_hi), 'n_T': int(t_n)},
        'phases_considered': phases,
        'phases_present_anywhere': present_any,
        'stability_ranges_K_by_phase': stability_ranges,
        'summary_by_T': summary_by_T,
        'ti2al3_like_phase_names': ti2al3_like,
        'endpoint_300K': {'stable_phases': stable_300K, 'fractions': fr_300K},
        'endpoint_1500K': {'stable_phases': stable_1500K, 'fractions': fr_1500K},
    }
    out_json.write_text(json.dumps(results, indent=2))

    # --- Print concise human summary ---
    print('Database:', tdb_path)
    print('Components:', components)
    print('Conditions: T[K]=[%.0f..%.0f] (%d pts), P=%.0f Pa, X(TI)=%.3f' % (t_lo, t_hi, t_n, p_pa, x_ti))
    print('Phases considered (%d): %s' % (len(phases), ', '.join(phases)))
    print('Phases present at any T:', ', '.join(present_any) if present_any else '(none)')
    print('Ti2Al3-like phase names (by pattern):', ', '.join(ti2al3_like) if ti2al3_like else '(none found)')
    print('Stable at 300 K:', ', '.join(f"{p}:{fr_300K[p]:.3f}" for p in sorted(fr_300K, key=lambda k: -fr_300K[k])) if fr_300K else '(none)')
    print('Stable at 1500 K:', ', '.join(f"{p}:{fr_1500K[p]:.3f}" for p in sorted(fr_1500K, key=lambda k: -fr_1500K[k])) if fr_1500K else '(none)')

    # Dominant phase and first-appearance temperature (≥10%)
    print('Dominant-phase appearance (≥10% fraction):')
    for pname in present_any:
        arr = np.asarray([float(frac_table[i, phase_names.index(pname)]) for i in range(len(Tvals))])
        idxs = np.where(arr >= 0.10)[0]
        if idxs.size:
            print(f"  {pname}: first ≥10% at ~{Tvals[idxs[0]]:.1f} K")
    print(f"Wrote: {out_json} and {out_csv}")


if __name__ == '__main__':
    main()
