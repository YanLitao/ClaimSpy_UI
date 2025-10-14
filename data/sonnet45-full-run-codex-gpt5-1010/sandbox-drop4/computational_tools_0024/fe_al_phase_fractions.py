#!/usr/bin/env python3
"""
Fe–Al CALPHAD: phase fractions near 50–52 at.% Al

Implements the single-script workflow per CALPHAD.md guidance:
- Load local TDB once (COST507R.TDB under TnETDBDB)
- Define conditions (T points and compositions)
- Compute equilibrium and phase fractions
- Summarize whether Al2Fe (FeAl2) or Fe2Al5 appear and quantify fractions
- Write results.json and results.csv; print concise human-readable summary

Notes
- Database path cached here from prior CALPHAD.md usage across this repo:
  /Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB
- Degrees of freedom: binary (AL, FE) → specify one v.X constraint
- Phase fractions computed from NP (phase moles) and normalized per condition
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from pycalphad import Database, equilibrium, variables as v


def phase_constituents(dbf: Database, phase: str) -> List[str]:
    const = dbf.phases[phase].constituents
    # Flatten unique constituents across all sublattices; convert Species→name
    flat = []
    for subl in const:
        for c in subl:
            name = getattr(c, 'name', str(c))
            flat.append(name)
    unique = sorted(set(flat))
    return unique


def select_phases(dbf: Database) -> List[str]:
    """
    Select phases that have any overlap with {AL, FE} in their constituents.
    This permissive filter avoids false negatives where TDB phases define
    multi-element sublattices but are still valid for AL–FE equilibria.
    """
    target = {"AL", "FE"}
    selected = []
    for ph in sorted(dbf.phases.keys()):
        try:
            const = set(phase_constituents(dbf, ph))
        except Exception:
            continue
        if const & target:
            selected.append(ph)
    return selected


def compute_phase_fractions_grouped(eq_ds):
    """
    Compute per-phase fractions at each temperature by grouping NP over 'vertex'
    using the 'Phase' labeling provided by pycalphad. Returns a list (over T) of
    dicts: [{phase_name: fraction, ...}, ...]. Assumes single composition and P.
    """
    pf = eq_ds["NP"].squeeze()          # dims: (T, vertex) after squeeze
    phase_lbl = eq_ds["Phase"].squeeze() # dims: (T, vertex)

    # Ensure 2D arrays [T, vertex]
    pf_np = np.asarray(pf.values)
    ph_np = np.asarray(phase_lbl.values)
    if pf_np.ndim == 1:
        pf_np = pf_np.reshape((1, -1))
    if ph_np.ndim == 1:
        ph_np = ph_np.reshape((1, -1))

    out: List[Dict[str, float]] = []
    for it in range(pf_np.shape[0]):
        phase_moles: Dict[str, float] = {}
        total = float(np.nansum(pf_np[it, :]))
        for iv in range(pf_np.shape[1]):
            name = str(ph_np[it, iv]) if ph_np[it, iv] is not None else ''
            name = name.strip()
            if not name:
                continue
            m = float(pf_np[it, iv])
            if not np.isfinite(m) or m <= 0:
                continue
            phase_moles[name] = phase_moles.get(name, 0.0) + m
        # Normalize to fractions
        if total > 0:
            for k in list(phase_moles.keys()):
                phase_moles[k] = phase_moles[k] / total
        out.append(phase_moles)
    return out


def first_non_nan_index(arr: np.ndarray) -> int:
    for i, val in enumerate(arr):
        if np.isfinite(val):
            return i
    return 0


def main():
    # Fixed inputs per prompt
    t_points = [298.0, 600.0, 800.0, 1000.0]
    x_al_list = [0.50, 0.51, 0.52]
    pressure = 101325.0  # Pa

    # Cache the TDB path (do not re-read CALPHAD.md repeatedly)
    tdb_path = Path("/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB")
    if not tdb_path.is_file():
        print(f"ERROR: TDB not found at {tdb_path}", file=sys.stderr)
        sys.exit(2)

    # Load database once
    dbf = Database(str(tdb_path))

    # Choose relevant phases
    phases = select_phases(dbf)
    if not phases:
        print("ERROR: No phases selected; check phase filter.", file=sys.stderr)
        sys.exit(3)

    # Target intermetallics (synonym-aware)
    key_targets = {
        "AL2FE": {"AL2FE", "FEAL2"},
        "AL5FE2": {"AL5FE2", "FE2AL5"},
    }

    # Compute and collect results
    results = {
        "database": str(tdb_path),
        "components": ["AL", "FE"],
        "phases_used": phases,
        "temperatures_K": t_points,
        "compositions_XAL": x_al_list,
        "entries": [],  # list of dicts per (X_AL, T)
    }

    csv_rows: List[Tuple[float, float, str, float]] = []

    for x_al in x_al_list:
        # Conditions: binary → one mole fraction specification
        conds = {
            v.T: t_points,
            v.P: pressure,
            v.N: 1.0,
            v.X("AL"): x_al,
        }

        try:
            eq_ds = equilibrium(dbf, ["AL", "FE"], phases, conds, output="NP")
        except Exception as e:
            print(f"ERROR during equilibrium at X(AL)={x_al}: {e}", file=sys.stderr)
            sys.exit(4)

        # Phase fractions via grouping by vertex labels
        temps = np.array(eq_ds["T"].values, dtype=float)
        grouped = compute_phase_fractions_grouped(eq_ds)

        # Build per-temperature entries
        for it, T in enumerate(temps):
            phase_fracs: Dict[str, float] = grouped[it]
            # Write CSV rows for present phases only
            for ph, val in phase_fracs.items():
                csv_rows.append((x_al, float(T), ph, float(val)))

            # Key targets
            targets_present = {}
            for label, synonyms in key_targets.items():
                present_val = 0.0
                for ph, val in phase_fracs.items():
                    if ph.upper() in synonyms:
                        present_val = max(present_val, val)
                targets_present[label] = present_val

            # Dominant phase
            dom_phase = max(phase_fracs.items(), key=lambda kv: kv[1])[0] if phase_fracs else ''

            results["entries"].append({
                "X_AL": x_al,
                "T": float(T),
                "dominant_phase": dom_phase,
                "targets": targets_present,
                "phase_fractions": phase_fracs,
            })

    # Write outputs
    out_json = Path("fe_al_results.json").resolve()
    out_csv = Path("fe_al_results.csv").resolve()

    with out_json.open("w") as f:
        json.dump(results, f, indent=2)

    with out_csv.open("w") as f:
        f.write("X_AL,T_K,phase,fraction\n")
        for x_al, T, ph, frac in csv_rows:
            f.write(f"{x_al:.4f},{T:.2f},{ph},{frac:.8f}\n")

    # Print concise human summary
    print(f"Using TDB: {tdb_path}")
    print(f"Phases considered ({len(phases)}): {', '.join(phases)}")
    print("")
    print("Summary (dominant phase and Al2Fe/Fe2Al5 fractions):")
    header = "X_AL    T[K]   Dominant         AL2FE_frac   AL5FE2_frac"
    print(header)
    print("-" * len(header))
    for entry in results["entries"]:
        x = entry["X_AL"]
        T = entry["T"]
        dom = entry["dominant_phase"]
        al2fe = entry["targets"].get("AL2FE", 0.0)
        al5fe2 = entry["targets"].get("AL5FE2", 0.0)
        print(f"{x:0.2f}  {T:6.1f}   {dom:12s}   {al2fe:11.4f}   {al5fe2:11.4f}")


if __name__ == "__main__":
    main()
