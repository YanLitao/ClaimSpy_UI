#!/usr/bin/env python3
"""
Equilibrium phase assemblage for Fe–Al at 750 K and 30 at.% Fe (70 at.% Al).

Workflow: load local TDB → select relevant phases → run pycalphad equilibrium →
compute phase fractions → write results.json and results.csv → print concise summary.

Database provenance: /Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB
Limitations per README: tuned for Al/Mg/Ti alloys but includes Fe binaries; suitable for Al–Fe.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf, phase, comps):
    const = dbf.phases[phase].constituents
    species = {c for subl in const for c in subl}
    return bool(set(comps) & species)


def main():
    # Fixed conditions (units: K, Pa, mole fraction)
    T = 750.0
    P = 101325.0
    x_al = 0.70  # 70 at.% Al → by balance 30 at.% Fe

    # Local database (absolute path as per CALPHAD.md guidance)
    tdb_path = "/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB"
    if not os.path.isfile(tdb_path):
        print(f"ERROR: TDB not found at {tdb_path}", file=sys.stderr)
        sys.exit(2)

    # Load DB once and introspect
    dbf = Database(tdb_path)

    # Components: include VA explicitly for substitutional solutions
    comps = ["AL", "FE", "VA"]

    # Quick sanity: ensure required elements exist in database
    # Prefer dbf.elements for element coverage
    elems = set(dbf.elements)
    missing = [c for c in ["AL", "FE"] if c not in elems]
    if missing:
        print(
            json.dumps(
                {
                    "status": "error",
                    "message": f"Database lacks required elements: {missing}",
                    "tdb": tdb_path,
                }
            )
        )
        sys.exit(1)

    # Phase list: include all phases in the database; single-point equilibrium is inexpensive
    phases = sorted(dbf.phases.keys())

    # Conditions: specify N−1 mole fractions for N substitutional components
    conds = {v.T: T, v.P: P, v.N: 1.0, v.X("AL"): x_al}
    # Mole-fraction constraints: exactly one is specified (binary system)

    # Compute equilibrium; use a modest pdens to aid convergence without excess cost
    eq = equilibrium(dbf, comps, phases, conds, calc_opts={"pdens": 1000})

    # DEBUG: print structure for robustness across pycalphad versions
    # (single run; acceptable lightweight introspection)
    try:
        print("EQ dims:", tuple(eq.dims))
        print("NP dims:", tuple(eq["NP"].dims))
        print("Phase dims:", tuple(eq["Phase"].dims))
    except Exception:
        pass

    # Phase moles → fractions by iterating vertices (robust across pycalphad versions)
    pf = eq["NP"]
    if "vertex" not in pf.dims:
        # Single-state case: treat total as one phase with unit fraction
        ph_val = eq["Phase"].values
        if hasattr(ph_val, "item"):
            ph_name = ph_val.item()
        else:
            ph_name = ph_val
        name = ph_name.decode() if isinstance(ph_name, (bytes, np.bytes_)) else str(ph_name)
        names = [name]
        ph_fracs = np.array([1.0], dtype=float)
    else:
        ph = eq["Phase"]
        nv = pf.sizes["vertex"]
        accum = {}
        total = 0.0
        for j in range(nv):
            np_j = float(np.squeeze(pf.isel(vertex=j).values))
            ph_j = ph.isel(vertex=j).values
            if hasattr(ph_j, "item"):
                ph_j = ph_j.item()
            name_j = ph_j.decode() if isinstance(ph_j, (bytes, np.bytes_)) else str(ph_j)
            if np_j > 0.0 and name_j not in ("", "<UNK>"):
                accum[name_j] = accum.get(name_j, 0.0) + np_j
                total += np_j
        # Normalize
        names = sorted(accum.keys(), key=lambda k: accum[k], reverse=True)
        if total > 0:
            ph_fracs = np.array([accum[k] / total for k in names], dtype=float)
        else:
            names = []
            ph_fracs = np.array([], dtype=float)

    # Build a dict of phase→fraction for those above a small threshold
    MINF = 1e-3
    results = {
        "T_K": float(T),
        "P_Pa": float(P),
        "composition": {"AL": float(x_al), "FE": float(1 - x_al)},
        "tdb_path": tdb_path,
        "phases": [],
    }
    for n, f in zip(names, ph_fracs):
        fv = float(f)
        if fv > MINF:
            results["phases"].append({"phase": n, "fraction": fv})

    # Sort by decreasing fraction for readability
    results["phases"].sort(key=lambda x: x["fraction"], reverse=True)

    # Write outputs to current directory
    out_json = Path("results.json")
    out_csv = Path("results.csv")
    with out_json.open("w") as f:
        json.dump(results, f, indent=2)

    with out_csv.open("w") as f:
        f.write("phase,fraction\n")
        for row in results["phases"]:
            f.write(f"{row['phase']},{row['fraction']:.6f}\n")

    # Human summary
    if not results["phases"]:
        print(
            "At 750 K and x_Al=0.70, no stable phases exceeded 1e-3 fraction (unexpected)."
        )
    else:
        print(
            f"Fe–Al at T=750 K, P=101325 Pa, x_Al=0.70 (x_Fe=0.30): {len(results['phases'])} phase(s)"
        )
        for row in results["phases"]:
            print(f"  {row['phase']}: {row['fraction']:.4f}")


if __name__ == "__main__":
    main()
