#!/usr/bin/env python3
"""
Fe–Al CALPHAD equilibrium at T=750 K, X_Fe=0.30 (X_Al=0.70).

Workflow:
- Load TDB from TnETDBDB directory
- Select phases that contain Fe/Al
- Compute equilibrium with pycalphad at fixed T,P,N and composition (N−1 mole-fraction constraints)
- Derive per‑phase mole fractions (from NP), list stable phases
- Verify the specific claim: two‑phase mixture of Al5Fe2 and Al2Fe at these conditions
- Write results.json and results.csv; print concise human summary

Notes:
- Degrees of freedom: for binary, set one v.X(...) constraint only
- Phase fractions computed as NP normalized by total moles; sum over vertex if present
"""

import json
import csv
from pathlib import Path
from typing import List, Dict

import numpy as np

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf: Database, phase: str, comps: List[str]) -> bool:
    const = dbf.phases[phase].constituents
    # Flatten constituents across sublattices; cast Species to str symbol names
    elems = {str(c).upper() for subl in const for c in subl}
    return bool(elems & {c.upper() for c in comps})


def main():
    cwd = Path.cwd()
    # TDB path (found once, not re-scanned repeatedly)
    tdb_path = Path("/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB")
    assert tdb_path.exists(), f"TDB not found: {tdb_path}"

    # System definition
    comps = ["AL", "FE"]
    T = 750.0  # K
    P = 101325.0  # Pa
    X_FE = 0.30
    # For binary, specify N−1 mole fraction constraints; X(AL) implied
    conds = {v.T: T, v.P: P, v.N: 1.0, v.X("FE"): X_FE}

    # Load database and select phases
    dbf = Database(str(tdb_path))
    phases = sorted(
        p for p in dbf.phases.keys() if phase_has_any(dbf, p, comps)
    )
    assert phases, "No phases selected; did you build the phases list?"

    # Compute equilibrium
    eq = equilibrium(dbf, comps, phases, conds)

    # Phase fractions (per-phase mole fraction vs T,P,composition)
    pf_da = eq["NP"]  # moles, per vertex if present
    fractions: Dict[str, float] = {}
    if "phase" in pf_da.dims:
        pf = pf_da
        if "vertex" in pf.dims:
            pf = pf.sum("vertex")
        total = float(pf.sum("phase").squeeze().values)
        phase_vals = pf.squeeze()
        phases_axis = phase_vals.coords["phase"].values.tolist()
        values = phase_vals.values
        if values.ndim != 1:
            other_axes = tuple(ax for ax in phase_vals.dims if ax != "phase")
            for ax in other_axes:
                phase_vals = phase_vals.squeeze(ax)
            phases_axis = phase_vals.coords["phase"].values.tolist()
            values = phase_vals.values
        for ph, val in zip(phases_axis, values):
            fractions[str(ph)] = float(val) / total if total > 0 else 0.0
    else:
        # Aggregate over vertices by phase name from the 'Phase' data variable
        # Squeeze singleton dims to length-vertex arrays
        pf = pf_da.squeeze()
        phase_names = eq["Phase"].squeeze()
        # Ensure we now have 1D arrays over 'vertex'
        if pf.ndim != 1 and "vertex" in pf.dims:
            other_axes = tuple(ax for ax in pf.dims if ax != "vertex")
            for ax in other_axes:
                pf = pf.squeeze(ax)
                phase_names = phase_names.squeeze(ax)
        # Map phase names to total moles
        totals: Dict[str, float] = {}
        for i in range(pf.sizes.get("vertex", pf.shape[0])):
            ph = str(phase_names.values[i])
            moles = float(pf.values[i])
            if ph:
                totals[ph] = totals.get(ph, 0.0) + moles
        total_all = sum(totals.values())
        fractions = {k: (v / total_all if total_all > 0 else 0.0) for k, v in totals.items()}

    # Filter to stable phases (fraction > threshold)
    thresh = 1e-3
    stable = {p: x for p, x in fractions.items() if x > thresh and np.isfinite(x)}

    # Verify the specific claim: two‑phase mixture of Al5Fe2 and Al2Fe
    def normalize(name: str) -> str:
        return name.strip().upper().replace(" ", "")

    targets = {"AL5FE2", "AL2FE"}
    stable_norm = {normalize(k): v for k, v in stable.items()}
    stable_set = set(stable_norm.keys())
    claim_ok = (len(stable_norm) == 2) and (targets <= stable_set)

    # Prepare outputs
    results = {
        "database": str(tdb_path),
        "components": comps,
        "conditions": {"T_K": T, "P_Pa": P, "X_FE": X_FE, "X_AL": 1 - X_FE},
        "phases_considered": phases,
        "stable_phases": [{"phase": k, "mole_fraction": stable[k]} for k in sorted(stable)],
        "all_phase_fractions": [{"phase": k, "mole_fraction": v} for k, v in sorted(fractions.items())],
        "claim": {
            "statement": "Equilibrium is a two-phase mixture of Al5Fe2 and Al2Fe",
            "valid": bool(claim_ok),
            "detail": {
                "stable_normalized": sorted(stable_set),
                "target_set": sorted(targets),
                "n_stable": len(stable_norm),
            },
        },
    }

    # Write JSON and CSV
    json_path = cwd / "fe_al_eq_results.json"
    with json_path.open("w") as fp:
        json.dump(results, fp, indent=2)

    csv_path = cwd / "fe_al_eq_results.csv"
    with csv_path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(["phase", "mole_fraction"])
        for row in sorted(stable.items()):
            writer.writerow(list(row))

    # Human summary
    print("Fe–Al equilibrium at T=750 K, X_Fe=0.30 (X_Al=0.70)")
    print(f"Database: {tdb_path}")
    if stable:
        print("Stable phases and mole fractions (threshold 1e-3):")
        for ph, mf in sorted(stable.items()):
            print(f"  - {ph}: {mf:.6f}")
    else:
        print("No stable phases above threshold.")
    print(f"Claim (two-phase Al5Fe2 + Al2Fe): {'VALID' if claim_ok else 'INVALID'}")
    print(f"Wrote: {json_path}")
    print(f"Wrote: {csv_path}")


if __name__ == "__main__":
    main()
