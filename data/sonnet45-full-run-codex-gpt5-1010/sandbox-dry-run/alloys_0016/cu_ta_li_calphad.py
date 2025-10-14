#!/usr/bin/env python3
"""
CALPHAD equilibrium for Cu–3Ta–0.5Li (at.%) from 298 K to 1073 K.

Implements a single-pass workflow:
- Load local COST507R.TDB from TnETDBDB
- Build relevant phase list by constituents
- Solve equilibrium vs T at fixed P and composition (N-1 X constraints)
- Derive per-phase fractions from NP and summarize stability
- Write results.json and results.csv and print concise human summary

Notes:
- Composition: Cu 0.965, Ta 0.030, Li 0.005 (atomic fractions)
- Temperature range: 298–1073 K (inclusive), 40 points (evenly spaced)
- Pressure: 101325 Pa; total moles N=1.0
- Phase fraction threshold for “present”: 1e-3
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf: Database, phase: str, comps: list[str]) -> bool:
    """Return True if the phase has any of the components in any sublattice."""
    const = dbf.phases[phase].constituents
    flat = {c for subl in const for c in subl}
    return bool(flat & set(comps))


def compute_equilibrium(tdb_path: Path,
                        temps: np.ndarray,
                        x_cu: float = 0.965,
                        x_ta: float = 0.030,
                        # x_li implied by balance
                        p_pa: float = 101325.0,
                        frac_thresh: float = 1e-3) -> dict:
    # Load database
    dbf = Database(str(tdb_path))

    # Basic checks: ensure elements exist
    elements = {el.upper() for el in dbf.elements}
    required = {"CU", "TA", "LI", "VA"}
    missing = required - elements
    if missing:
        raise RuntimeError(f"Selected TDB lacks required elements: {sorted(missing)}")

    # Components and phase selection (constituent filter)
    comps = ["CU", "TA", "LI", "VA"]

    # Filter phases to those that include any of the active components (excluding VA-only)
    phases = sorted(p for p in dbf.phases if phase_has_any(dbf, p, ["CU", "TA", "LI"]))
    if not phases:
        # Fallback: consider all phases; some databases list constituents per phase narrowly
        phases = sorted(dbf.phases.keys())

    # Conditions: N-1 mole fraction constraints
    if not np.isclose(x_cu + x_ta, 0.995):
        raise ValueError("Composition mismatch: X(CU)+X(TA) must be 0.995 for 0.5 at.% Li.")

    conds = {
        v.T: temps,
        v.P: p_pa,
        v.N: 1.0,
        v.X("CU"): x_cu,
        v.X("TA"): x_ta,
        # X(LI) implied by balance in a 3-component substitutional set
    }

    # Solve equilibrium
    eq = equilibrium(dbf, comps, phases, conds, verbose=False)
    if "NP" not in eq.data_vars:
        raise RuntimeError(f"Equilibrium result lacks 'NP'; available vars: {list(eq.data_vars)}")

    # Compute per-phase fractions
    pf = eq["NP"]  # moles of each phase per condition
    # In recent pycalphad, 'phase' is typically a coordinate on 'vertex' dimension; groupby works
    # Build per-phase moles while preserving condition dims by summing over 'vertex'
    if "vertex" in pf.dims:
        # Prefer a label variable where available
        label_var = None
        for cand in ["phase", "Phase", "PHASE", "phase_name", "PHASE_NAME"]:
            if cand in eq and "vertex" in eq[cand].dims:
                label_var = cand
                break
        if label_var is None:
            raise RuntimeError(
                f"Cannot find a phase-label variable; dims={pf.dims}, coords={list(pf.coords)}, vars={list(eq.data_vars)}"
            )
        label = eq[label_var]
        unique_phases = sorted(str(x) for x in np.unique(label.values))
        pieces = []
        for ph in unique_phases:
            sel = pf.where(label == ph)
            arr = sel.fillna(0).sum("vertex")
            arr = arr.expand_dims({"phase": [ph]})
            pieces.append(arr)
        pf_by_phase = xr.concat(pieces, dim="phase")
    elif "phase" in pf.dims:
        pf_by_phase = pf
    else:
        raise RuntimeError(
            f"Cannot interpret 'NP' dimensions; dims={pf.dims}, coords={list(pf.coords)}, vars={list(eq.data_vars)}"
        )
    # Normalize per-phase fractions
    tot = pf_by_phase.sum("phase")
    phase_frac = (pf_by_phase / tot)

    # Reduce to T-major table for CSV/JSON
    # The dataset dims typically include (T, P, X_CU, X_TA, N, phase)
    # We collapse all but T and phase since they are scalars here.
    # Convert to tidy DataFrame then pivot to wide (T x phase)
    df_long = phase_frac.to_dataframe(name="fraction").reset_index()
    # Keep only T and phase; average over any degenerate dims if present
    if "T" not in df_long.columns or "phase" not in df_long.columns:
        raise RuntimeError(f"Unexpected columns in long DF: {df_long.columns.tolist()}")
    df_grouped = df_long.groupby(["T", "phase"], as_index=False)["fraction"].mean()
    df = df_grouped.pivot(index="T", columns="phase", values="fraction").fillna(0.0)
    df.index.name = "T_K"
    df.sort_index(inplace=True)

    # Determine stable phases (threshold)
    stable_mask = (df > frac_thresh)

    # Endpoint phases
    t_min, t_max = float(df.index.min()), float(df.index.max())
    phases_at_min = sorted(df.columns[stable_mask.loc[t_min]]) if t_min in stable_mask.index else []
    phases_at_max = sorted(df.columns[stable_mask.loc[t_max]]) if t_max in stable_mask.index else []

    # Dominant matrix phase at high-T (by fraction at t_max)
    top_phase = None
    if t_max in df.index:
        top_phase = df.loc[t_max].idxmax()

    # First-appearance temperature for non-matrix phases
    precip_info: dict[str, float] = {}
    for phase in df.columns:
        if phase == top_phase:
            continue
        present = stable_mask[phase]
        if present.any():
            # first T where present is True
            first_idx = present.idxmax() if present.iloc[0] else present.idxmax()
            # idxmax returns first occurrence of max(True==1); guard by checking any
            t_first = float(first_idx)
            precip_info[phase] = t_first

    # Single-phase vs multi-phase determination across the range
    num_phases_over_T = stable_mask.sum(axis=1)
    single_phase_all = bool((num_phases_over_T <= 1).all())

    # Summaries
    summary = {
        "tdb": str(tdb_path),
        "components": comps,
        "phases_considered": phases,
        "temperature_K_range": [t_min, t_max],
        "composition_at_frac": {"X(CU)": x_cu, "X(TA)": x_ta, "X(LI)": 1.0 - x_cu - x_ta},
        "frac_threshold": frac_thresh,
        "stable_phases_at_298K": phases_at_min,
        "stable_phases_at_1073K": phases_at_max,
        "matrix_phase_highT": top_phase,
        "precipitates_first_appearance_K": precip_info,
        "single_phase_all_T": single_phase_all,
    }

    return {
        "eq_table": df,
        "summary": summary,
    }


def main() -> int:
    # Inputs
    tdb_path = Path("/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB").resolve()
    temps = np.linspace(298.0, 1073.0, 40)

    try:
        results = compute_equilibrium(tdb_path=tdb_path, temps=temps)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    df = results["eq_table"]
    summary = results["summary"]

    # Write outputs
    out_json = Path("results.json")
    out_csv = Path("results.csv")
    df.to_csv(out_csv)
    with open(out_json, "w") as f:
        json.dump({
            "summary": summary,
            "temperatures_K": df.index.tolist(),
            "phase_fractions": df.to_dict(orient="list"),
        }, f, indent=2)

    # Human summary
    comp = summary["composition_at_frac"]
    print("CALPHAD summary — Cu–3Ta–0.5Li (at.%)")
    print(f"TDB: {summary['tdb']}")
    print(f"T-range: {summary['temperature_K_range'][0]:.1f}–{summary['temperature_K_range'][1]:.1f} K; P={101325} Pa")
    print(f"Composition: X(CU)={comp['X(CU)']:.4f}, X(TA)={comp['X(TA)']:.4f}, X(LI)={comp['X(LI)']:.4f}")
    print(f"Phases considered (n={len(summary['phases_considered'])}): {', '.join(summary['phases_considered'])}")
    print(f"Stable at 298 K: {', '.join(summary['stable_phases_at_298K']) if summary['stable_phases_at_298K'] else 'None above threshold'}")
    print(f"Stable at 1073 K: {', '.join(summary['stable_phases_at_1073K']) if summary['stable_phases_at_1073K'] else 'None above threshold'}")
    print(f"High-T matrix phase: {summary['matrix_phase_highT']}")
    if summary['precipitates_first_appearance_K']:
        pairs = [f"{k} @ {v:.1f} K" for k, v in sorted(summary['precipitates_first_appearance_K'].items(), key=lambda kv: kv[1])]
        print("Precipitates (first appearance): " + "; ".join(pairs))
    else:
        print("Precipitates (first appearance): none above threshold")
    print(f"Single-phase across range: {summary['single_phase_all_T']}")
    print(f"Wrote: {out_csv} and {out_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
