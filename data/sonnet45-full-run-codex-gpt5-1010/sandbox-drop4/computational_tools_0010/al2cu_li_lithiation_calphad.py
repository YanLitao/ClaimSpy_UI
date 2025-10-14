#!/usr/bin/env python3
"""
CALPHAD study: Equilibrium phase assemblage during initial lithiation of Al2Cu.

Assumptions and setup
- Composition path: start from stoichiometric Al2Cu (Al:Cu = 2:1) and add Li from 0 to 0.15 at. fraction of the total.
  Thus X_LI = x in [0, 0.15]; X_AL = (1 − x) * 2/3; X_CU = (1 − x) * 1/3.
- Temperatures: 298 K and 400 K; pressure 1 atm (101325 Pa); total moles normalized to 1.
- Database: COST507R (local), loaded from TnETDBDB directory.
- Computation: sequential equilibrium solves for each composition (avoid Cartesian product meshes).
- Reporting: stable phases and phase fractions at each (T, x_Li). Initial lithiation region defined here as 0.00–0.10 Li.

Outputs
- results JSON: al2cu_li_lithiation_results.json
- results CSV:  al2cu_li_lithiation_results.csv
- Printed human summary including initial lithiation classification (single-/two-/multi-phase) at each T.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pycalphad import Database, equilibrium, variables as v


P_ATM = 101325
NP_TOL = 1e-6


def find_tdb_candidates(tdb_dir: str) -> List[str]:
    """Return absolute paths to .tdb/.TDB files in tdb_dir (non-recursive)."""
    p = Path(tdb_dir)
    cands = []
    for child in p.iterdir():
        if child.is_file() and child.suffix.lower() == ".tdb":
            cands.append(str(child.resolve()))
    return sorted(cands)


def select_db_for_elements(candidates: List[str], required: List[str]) -> Tuple[str, List[str]]:
    """Load each candidate once and select one whose elements superset 'required'. Return chosen path and considered list.
    Preference heuristic: keep the first that covers required elements; tie-break by filename containing 'cost507' if present.
    """
    considered = []
    viable = []
    for path in candidates:
        considered.append(path)
        try:
            db = Database(path)
            elems = set(db.elements)
            if set(required).issubset(elems):
                viable.append((path, elems))
        except Exception:
            continue
    if not viable:
        return "", considered
    # Prefer COST 507 if present
    viable_sorted = sorted(viable, key=lambda t: ("cost507" not in os.path.basename(t[0]).lower(), t[0]))
    return viable_sorted[0][0], considered


def phases_allowed_only(db: Database, phase: str, allowed_species: List[str]) -> bool:
    """True if all constituents of 'phase' are subset of allowed_species."""
    allowed = set(allowed_species)
    pobj = db.phases[phase]
    species = set()
    for subl in pobj.constituents:
        species.update(subl)
    return species.issubset(allowed)


def curated_phase_lists(db: Database) -> Tuple[List[str], List[str]]:
    """Curate two phase lists from database names:
    - Binary Al–Cu list (no Li-bearing compounds)
    - Ternary Al–Cu–Li list (includes Li-bearing intermetallics)
    Always include common solution phases and LIQUID if present.
    """
    common = [p for p in ["LIQUID", "FCC_A1", "BCC_A2", "HCP_A3"] if p in db.phases]
    alcu = sorted([p for p in db.phases if p.startswith("ALCU_")])
    alculi = sorted([p for p in db.phases if p.startswith("ALCULI_")])
    # Binary excludes Li-bearing intermetallics
    phases_bin = sorted(set(common + alcu))
    # Ternary includes both
    phases_ternary = sorted(set(common + alcu + alculi))
    if not phases_bin:
        raise RuntimeError("No suitable binary Al–Cu phases found in DB.")
    if not phases_ternary:
        raise RuntimeError("No suitable ternary Al–Cu–Li phases found in DB.")
    return phases_bin, phases_ternary


def compute_equilibrium_series(db: Database,
                               phases: List[str],
                               temps: List[float],
                               x_li_grid: np.ndarray,
                               out_json: str,
                               out_csv: str,
                               tdb_path: str) -> Dict:
    comps = ["AL", "CU", "LI", "VA"]
    results_rows = []
    summary = {
        "tdb_path": tdb_path,
        "components": comps,
        "phases": phases,
        "pressure_Pa": P_ATM,
        "np_tol": NP_TOL,
        "temps_K": temps,
        "x_li_grid": x_li_grid.tolist(),
        "points": []  # list of dicts per (T, x_li)
    }

    for T in temps:
        for x_li in x_li_grid:
            # Composition mapping from Al2Cu baseline
            x_al = (1.0 - x_li) * (2.0 / 3.0)
            x_cu = (1.0 - x_li) * (1.0 / 3.0)

            conds = {v.T: float(T), v.P: P_ATM, v.N: 1.0, v.X("AL"): float(x_al), v.X("CU"): float(x_cu)}

            # Guard degrees of freedom: ternary → 2 mole-fraction constraints provided
            if abs((x_al + x_cu + x_li) - 1.0) > 1e-8:
                raise ValueError("Composition mapping error: fractions do not sum to unity.")

            try:
                # Use binary Al–Cu phase list when x_li is effectively zero to avoid degeneracy
                if x_li <= 0.0 + 1e-12:
                    phases_bin, phases_ternary = curated_phase_lists(db)
                    # Compute as a binary subsystem by excluding Li from components
                    eq = equilibrium(db, ["AL", "CU", "VA"], phases_bin,
                                     {v.T: float(T), v.P: P_ATM, v.N: 1.0, v.X("AL"): float(2/3), v.X("CU"): float(1/3)},
                                     calc_opts={"pdens": 1200})
                else:
                    # Use curated ternary list
                    _, phases_ternary = curated_phase_lists(db)
                    eq = equilibrium(db, comps, phases_ternary, conds, calc_opts={"pdens": 1200})
            except Exception:
                # Retry with solution phases only
                simple = [p for p in ["LIQUID", "FCC_A1", "BCC_A2", "HCP_A3"] if p in db.phases]
                eq = equilibrium(db, comps, simple, conds)

            # Phase moles are attached to vertices; aggregate by Phase label per-vertex
            # Shapes: NP dims ~ (N,P,T,X_AL,X_CU,vertex); Phase dims ~ (N,P,T,X_AL,X_CU,vertex)
            np_da = eq["NP"].values  # ndarray
            ph_da = eq["Phase"].values  # ndarray of dtype object/bytes
            # Squeeze condition dims -> 1D over vertex
            np_vec = np.squeeze(np_da)
            ph_vec = np.squeeze(ph_da)
            # Ensure 1D vertex vectors
            if np_vec.ndim == 0:
                np_vec = np.array([float(np_vec)])
                ph_vec = np.array([ph_vec])
            assert np_vec.ndim == 1 and ph_vec.ndim == 1 and np_vec.shape[0] == ph_vec.shape[0]

            # Group by phase name and sum moles across vertices
            agg: Dict[str, float] = {}
            for a, p in zip(np_vec, ph_vec):
                name = p.decode() if isinstance(p, (bytes, bytearray)) else str(p)
                aval = float(a)
                if not np.isfinite(aval) or aval <= 0.0:
                    continue
                agg[name] = agg.get(name, 0.0) + aval

            # Normalize by total to get fractions
            tot = sum(agg.values())
            present = []
            if tot > 0:
                present = sorted(((k, v / tot) for k, v in agg.items() if v / tot > 1e-6), key=lambda t: -t[1])

            results_rows.append({
                "T_K": T,
                "x_LI": float(x_li),
                "x_AL": float(x_al),
                "x_CU": float(x_cu),
                "num_phases": len(present),
                "phases": [p for p, _ in present],
                "fractions": [f for _, f in present]
            })

            summary["points"].append({
                "T_K": T,
                "x_LI": float(x_li),
                "phases": [{"name": p, "fraction": f} for p, f in present]
            })

    # Write CSV
    df = pd.DataFrame(results_rows)
    # For compact CSV, add a semicolon-delimited phase:frac column
    def fmt_pairs(row):
        return "; ".join(f"{p}:{f:.6f}" for p, f in zip(row["phases"], row["fractions"]))

    df["assemblage"] = df.apply(fmt_pairs, axis=1)
    df_out = df[["T_K", "x_LI", "x_AL", "x_CU", "num_phases", "assemblage"]].copy()
    df_out.to_csv(out_csv, index=False)

    # Write JSON
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def classify_initial_region(points: List[Dict], T: float, x_lo: float = 0.0, x_hi: float = 0.10,
                            minf: float = 1e-3) -> Dict:
    """For a given T, inspect points within x_lo..x_hi and classify the assemblage multiplicity.
    Returns a dict with counts of unique phase-counts observed and a suggested label.
    """
    pts = [p for p in points if abs(p["T_K"] - T) < 1e-9 and x_lo - 1e-12 <= p["x_LI"] <= x_hi + 1e-12]
    if not pts:
        return {"T_K": T, "classification": "insufficient-data", "details": {}}

    counts = []
    examples = {}
    for p in pts:
        # count phases above threshold
        n = sum(1 for ph in p["phases"] if ph["fraction"] > minf)
        counts.append(n)
        examples.setdefault(n, []).append((p["x_LI"], [(ph["name"], ph["fraction"]) for ph in p["phases"] if ph["fraction"] > minf]))

    unique_counts = sorted(set(counts))
    if len(unique_counts) == 1:
        n = unique_counts[0]
        label = {1: "single-phase", 2: "two-phase"}.get(n, "multi-phase")
    else:
        # mixture across the region
        if max(unique_counts) <= 2:
            label = "mixed single-/two-phase"
        else:
            label = "mixed multi-phase"

    return {"T_K": T, "classification": label, "phase_count_values": unique_counts, "examples": examples}


def main():
    # Locate TDB directory from project root (absolute), per CALPHAD.md
    TDB_DIR = "/Users/delip/play/miniclaimspy/TnETDBDB"
    required_elems = ["AL", "CU", "LI", "VA"]

    cands = find_tdb_candidates(TDB_DIR)
    chosen, considered = select_db_for_elements(cands, required_elems)
    if not chosen:
        print("ERROR: No local TDB contains the required elements AL, CU, LI.")
        print("Considered:")
        for c in considered:
            print(f"  - {c}")
        raise SystemExit(2)

    # Load chosen DB and shortlist phases
    db = Database(chosen)
    # Curate phase lists
    phases_bin, phases_ternary = curated_phase_lists(db)
    phases = phases_ternary

    # Composition grid: Li from 0 to 0.15 (inclusive) at fine steps
    x_li_grid = np.round(np.linspace(0.0, 0.15, 16), 4)
    temps = [298.0, 400.0]

    out_json = os.path.join(os.getcwd(), "al2cu_li_lithiation_results.json")
    out_csv = os.path.join(os.getcwd(), "al2cu_li_lithiation_results.csv")

    summary = compute_equilibrium_series(db, phases, temps, x_li_grid, out_json, out_csv, chosen)

    # Classify initial lithiation region (0–10 at.% Li)
    print("Database:", summary["tdb_path"]) 
    print("Temperatures:", ", ".join(f"{T:.1f} K" for T in temps))
    print("Li grid (at. fraction):", ", ".join(f"{x:.3f}" for x in x_li_grid))
    print()

    for T in temps:
        cls = classify_initial_region(summary["points"], T, 0.01, 0.10, minf=1e-3)
        print(f"Initial lithiation (0–10 at.% Li) at {T:.0f} K: {cls['classification']}")
        # Show a few representative compositions
        # pick x in {0.00, 0.02, 0.05, 0.10}
        targets = [0.00, 0.02, 0.05, 0.10]
        for xt in targets:
            # find closest available
            pt = min((p for p in summary["points"] if abs(p["T_K"] - T) < 1e-9), key=lambda p: abs(p["x_LI"] - xt))
            assemblage = ", ".join(f"{ph['name']}:{ph['fraction']:.3f}" for ph in pt["phases"] if ph["fraction"] > 1e-3)
            print(f"  x(Li)≈{pt['x_LI']:.3f}: {assemblage}")
        print()

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_json}")


if __name__ == "__main__":
    main()
