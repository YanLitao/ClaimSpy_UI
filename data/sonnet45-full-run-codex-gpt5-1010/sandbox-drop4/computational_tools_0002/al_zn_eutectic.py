#!/usr/bin/env python3
"""
CALPHAD analysis for the Al–Zn binary to identify the eutectic point.

Follows workflow constraints:
- Uses only local TDBs from TnETDBDB (repo-provided DBs).
- Single script: load DB → set conditions → solve → summarize → write JSON/CSV.
- Fix DOF: N−1 mole-fraction constraints for a binary.
- Summarize stable phases and identify eutectic via global minimum of the liquidus.

Outputs:
- results.json: structured summary (eutectic T, composition in at% Zn, DB used)
- results.csv: liquidus T vs composition (coarse and refined if available)
- liquidus.png: visualization saved (no interactive show)
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from dataclasses import dataclass
from typing import List, Dict, Any

import numpy as np
import pandas as pd

from pycalphad import Database, variables as v, equilibrium


# Cache the path hints from CALPHAD.md without rereading during runtime.
# Repository convention: local TDBs live under this directory.
TDB_DIR = "/Users/delip/play/miniclaimspy/TnETDBDB"


@dataclass
class EutecticResult:
    tdb_path: str
    components: List[str]
    phases_considered: List[str]
    eutectic_T_K: float | None
    eutectic_x_zn: float | None
    eutectic_atpct_zn: float | None
    message: str


def list_tdb_files(tdb_dir: str) -> List[str]:
    if not os.path.isdir(tdb_dir):
        return []
    files = []
    for fn in os.listdir(tdb_dir):
        if fn.lower().endswith('.tdb'):
            files.append(os.path.join(tdb_dir, fn))
    return sorted(files)


def phase_has_any(dbf: Database, phase: str, comps: List[str]) -> bool:
    const = dbf.phases[phase].constituents
    present = {str(c).upper() for subl in const for c in subl}
    return bool(present & set(comps))


def select_db_for_al_zn(tdb_candidates: List[str]) -> str | None:
    """Select the first database that contains both AL and ZN.
    Prefer commonly used Al-alloy DBs if present (e.g., COST507*).
    """
    # Preferential ordering: any file starting with cost507 (case-insensitive)
    preferred = [p for p in tdb_candidates if os.path.basename(p).lower().startswith('cost507')]
    ordered = preferred + [p for p in tdb_candidates if p not in preferred]
    for path in ordered:
        try:
            dbf = Database(path)
            syms = {str(s).upper() for s in dbf.elements}
            if {'AL', 'ZN', 'VA'}.issubset(syms):
                return path
        except Exception:
            continue
    return None


def compute_equilibrium(dbf: Database, phases: List[str], T: np.ndarray, X_zn: np.ndarray):
    comps = ['AL', 'ZN', 'VA']
    # Conditions per CALPHAD workflow: N−1 composition constraints for N substitutional comps
    conds = {
        v.T: T,
        v.P: 101325.0,
        v.N: 1.0,
        v.X('ZN'): X_zn,  # X(AL) is implied
    }
    # Compute equilibrium; suppress warnings to keep logs clean
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        eq = equilibrium(dbf, comps, phases, conds, verbose=False)
    return eq


def phase_fractions(eq) -> Any:
    """Return NP summed over 'vertex' and normalized by total at each (N,P,T,X).
    Note: Use eq['Phase'] alongside this to identify which entries correspond to which phase.
    """
    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')
    # Total moles per (N,P,T,X)
    tot = pf.sum(dim=('N', 'P'))
    return pf / tot


def liquidus_min_T(eq, thresh=1e-8):
    """Minimal temperature where LIQUID participates in equilibrium vs composition.

    Uses eq['Phase'] (phase labels over N,P,vertex) and eq['NP'] (moles) to detect
    whether LIQUID is present above a small threshold at each (T, X_ZN).
    Returns (X_ZN array, T_liquidus array with NaN where no liquid in window).
    """
    T_vals = np.array(eq['T'].values)
    X_vals = np.array(eq['X_ZN'].values)
    np_arr = np.array(eq['NP'].values)        # (N, P, T, X, vertex)
    ph_arr = np.array(eq['Phase'].values)     # same shape
    assert np_arr.shape == ph_arr.shape
    is_liq = (ph_arr == 'LIQUID')
    has_liq = (np_arr > thresh) & is_liq
    presence = has_liq.any(axis=(0, 1, 4))    # (T, X)
    nT, nX = presence.shape
    T_min = np.full(nX, np.nan, dtype=float)
    for j in range(nX):
        if presence[:, j].any():
            first_true = np.where(presence[:, j])[0][0]
            T_min[j] = T_vals[first_true]
    return X_vals, T_min


def summarize_eutectic(T_vals: np.ndarray, X_vals: np.ndarray, T_liq_min: np.ndarray) -> Dict[str, Any]:
    # Identify global minimum of liquidus temperature across composition
    finite = np.isfinite(T_liq_min)
    if not finite.any():
        return {
            'eutectic_T_K': None,
            'eutectic_x_zn': None,
            'eutectic_atpct_zn': None,
            'message': 'Liquidus not found in sampled window.'
        }
    jstar = int(np.nanargmin(T_liq_min))
    T_star = float(T_liq_min[jstar])
    x_star = float(X_vals[jstar])
    return {
        'eutectic_T_K': T_star,
        'eutectic_x_zn': x_star,
        'eutectic_atpct_zn': 100.0 * x_star,
        'message': 'Global minimum of liquidus across composition.'
    }


def main() -> int:
    # Discover TDBs once
    tdb_files = list_tdb_files(TDB_DIR)
    if not tdb_files:
        res = EutecticResult(
            tdb_path='', components=['AL', 'ZN', 'VA'], phases_considered=[],
            eutectic_T_K=None, eutectic_x_zn=None, eutectic_atpct_zn=None,
            message='No .TDB files found in TnETDBDB; cannot perform Al–Zn analysis.'
        )
        with open('results.json', 'w') as f:
            json.dump(res.__dict__, f, indent=2)
        print(res.message)
        return 0

    tdb_path = select_db_for_al_zn(tdb_files)
    if tdb_path is None:
        res = EutecticResult(
            tdb_path='', components=['AL', 'ZN', 'VA'], phases_considered=[],
            eutectic_T_K=None, eutectic_x_zn=None, eutectic_atpct_zn=None,
            message='Available local TDBs do not provide required elements {AL, ZN, VA}.'
        )
        with open('results.json', 'w') as f:
            json.dump(res.__dict__, f, indent=2)
        print(res.message)
        return 0

    # Load DB once
    dbf = Database(tdb_path)
    # Prefer a focused set of phases for Al–Zn
    candidate_core = ['LIQUID', 'FCC_A1', 'HCP_A3', 'HCP_ZN']
    phases = [p for p in candidate_core if p in dbf.phases]
    # Fallback: all phases containing AL or ZN
    if not phases:
        phases = sorted(p for p in dbf.phases if phase_has_any(dbf, p, ['AL', 'ZN']))
    assert phases, "No phases selected; did you build the phases list?"

    # Debug: report LIQUID constituents once (useful to understand availability)
    try:
        const = dbf.phases['LIQUID'].constituents
        liq_const = sorted({str(c).upper() for subl in const for c in subl})
        print("LIQUID constituents:", ', '.join(liq_const))
    except Exception:
        print("LIQUID phase not present in database constituents (unexpected).")

    # Coarse sweep
    T_coarse = np.linspace(450.0, 950.0, 81)
    X_coarse = np.linspace(0.0, 1.0, 121)
    eq_c = compute_equilibrium(dbf, phases, T_coarse, X_coarse)
    # Minimal debug: list variables and their dims once
    try:
        var_info = {k: tuple(v.dims) for k, v in eq_c.data_vars.items()}
        print("DEBUG eq vars:", {k: var_info[k] for k in sorted(var_info) if k in ['NP','GM','MU','X','Y']})
        print("DEBUG eq coords:", list(eq_c.coords))
    except Exception:
        pass
    pf_c = phase_fractions(eq_c)
    # Debug: coarse liquid presence check via Phase/NP
    try:
        Xc_dbg, Tliq_dbg = liquidus_min_T(eq_c, thresh=1e-8)
        print("DEBUG liquidus sample:",
              "found" if np.isfinite(np.nanmin(Tliq_dbg)) else "none in window")
    except Exception:
        pass
    Xc, Tliq_min_c = liquidus_min_T(eq_c, thresh=1e-8)
    coarse_summary = summarize_eutectic(T_coarse, Xc, Tliq_min_c)

    # Optionally refine near the coarse minimum
    T_star = coarse_summary['eutectic_T_K']
    x_star = coarse_summary['eutectic_x_zn']
    if T_star is not None and x_star is not None:
        dT = 60.0
        dX = 0.12
        T_ref = np.linspace(max(350.0, T_star - dT), min(1000.0, T_star + dT), 121)
        X_ref = np.linspace(max(0.0, x_star - dX), min(1.0, x_star + dX), 241)
        eq_r = compute_equilibrium(dbf, phases, T_ref, X_ref)
        pf_r = phase_fractions(eq_r)
        Xr, Tliq_min_r = liquidus_min_T(eq_r, thresh=1e-8)
        refined_summary = summarize_eutectic(T_ref, Xr, Tliq_min_r)
        # Prefer refined if it produced a lower liquidus minimum
        if (refined_summary['eutectic_T_K'] is not None and
            refined_summary['eutectic_T_K'] <= coarse_summary['eutectic_T_K'] + 1e-6):
            best = refined_summary
            best_T_grid = T_ref
            best_X_grid = Xr
            best_Tliquidus = Tliq_min_r
        else:
            best = coarse_summary
            best_T_grid = T_coarse
            best_X_grid = Xc
            best_Tliquidus = Tliq_min_c
    else:
        best = coarse_summary
        best_T_grid = T_coarse
        best_X_grid = Xc
        best_Tliquidus = Tliq_min_c

    # Prepare outputs
    res = EutecticResult(
        tdb_path=tdb_path,
        components=['AL', 'ZN', 'VA'],
        phases_considered=phases,
        eutectic_T_K=best['eutectic_T_K'],
        eutectic_x_zn=best['eutectic_x_zn'],
        eutectic_atpct_zn=best['eutectic_atpct_zn'],
        message=best['message'],
    )

    # Verify phase assemblages just below/above eutectic (if found)
    phase_check = {}
    if res.eutectic_T_K is not None and res.eutectic_x_zn is not None:
        T_probe = [max(1.0, res.eutectic_T_K - 2.0), res.eutectic_T_K + 2.0]
        X_probe = [min(0.999999, max(1e-9, res.eutectic_x_zn))]
        eq_probe = compute_equilibrium(dbf, phases, np.array(T_probe), np.array(X_probe))
        np_arr = np.array(eq_probe['NP'].values)
        ph_arr = np.array(eq_probe['Phase'].values)
        mask = np_arr > 1e-8
        # Below eutectic (index 0) and above eutectic (index 1)
        for i, label in enumerate(['below', 'above']):
            present = ph_arr[:, :, i, 0, :][mask[:, :, i, 0, :]]
            if present.size == 0:
                phase_check[label] = []
            else:
                phase_check[label] = sorted(set(str(p) for p in present))
    

    # Write JSON
    # Write JSON (include phase_check if computed)
    out = res.__dict__.copy()
    if 'phase_check' in locals():
        out['phase_check'] = phase_check
    with open('results.json', 'w') as f:
        json.dump(out, f, indent=2)

    # Write CSV (liquidus curve)
    df = pd.DataFrame({
        'X_ZN': best_X_grid,
        'T_liquidus_min_K': best_Tliquidus,
        'atpct_Zn': 100.0 * best_X_grid
    })
    df.to_csv('results.csv', index=False)

    # Plot liquidus as a check (optional)
    try:
        import matplotlib
        matplotlib.use('Agg')  # no display
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 4))
        plt.plot(100.0 * best_X_grid, best_Tliquidus, '-', lw=1.5, color='tab:blue', label='Liquidus (min T vs X)')
        if res.eutectic_atpct_zn is not None and res.eutectic_T_K is not None:
            plt.scatter([res.eutectic_atpct_zn], [res.eutectic_T_K], color='red', zorder=5, label='Eutectic (min)')
        plt.xlabel('Zn (at%)')
        plt.ylabel('Temperature (K)')
        plt.title('Al–Zn Liquidus (from local TDB)')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig('liquidus.png', dpi=160)
    except Exception:
        pass

    # Concise human summary
    print("Database:", res.tdb_path)
    print("Phases considered:", ', '.join(res.phases_considered))
    if res.eutectic_T_K is not None and res.eutectic_atpct_zn is not None:
        print(f"Eutectic (by liquidus minimum): T = {res.eutectic_T_K:.1f} K, Zn = {res.eutectic_atpct_zn:.2f} at%")
        if 'phase_check' in locals():
            print("Assemblages near eutectic:", phase_check)
    else:
        print("Eutectic not found in the sampled window.")

    return 0


if __name__ == '__main__':
    sys.exit(main())
