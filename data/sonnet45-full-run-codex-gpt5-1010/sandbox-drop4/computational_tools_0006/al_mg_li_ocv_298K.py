"""
Compute equilibrium open-circuit potential (OCV) vs Li/Li+ for Li insertion
into an Al–Mg alloy anode at 298 K using CALPHAD (pycalphad).

Workflow (per cached essentials from CALPHAD.md):
- Use only local TDBs in TnETDBDB; here COST507R.TDB at
  /Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB
- Single script: load DB → set conditions → solve → summarize → write JSON/CSV.
- Fix degrees of freedom and fractions: v.N = 1.0; specify N−1 mole fractions.
- Use pycalphad to list phases and perform equilibrium; avoid repeated DB loads.

Thermodynamic formulation for OCV vs Li/Li+:
  E = - (mu_Li^(alloy) - mu_Li^(metal)) / F
where mu_Li^(alloy) is the equilibrium chemical potential of Li in the
Al–Mg–Li system at the chosen bulk composition and temperature, and
mu_Li^(metal) is the chemical potential of Li in pure Li metal at 298 K.
F = 96485 C/mol.

We evaluate mu_Li^(alloy) at 298 K for a series of small Li contents (x_Li)
while preserving the Al:Mg ratio at 1:1 (host ≈ Al0.5Mg0.5 before insertion):
  X(LI) = x
  X(AL) = 0.5 * (1 - x)
  X(MG) is implied from the constraint that mole fractions sum to 1.

Outputs:
- results.json: numerical results and metadata
- summary.csv: table of x_Li, mu_Li(alloy), mu_Li(metal), OCV(V), dominant phases
- summary.txt: concise human-readable summary

Note: The COST507 Al database includes Li, Al, Mg and their intermetallics
(e.g., ALLI, AL2LI3, AL12MG17) and solution phases (FCC_A1, BCC_A2, HCP_A3).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Dict, Any

import numpy as np

from pycalphad import Database, equilibrium, variables as v


FARADAY = 96485.0  # C/mol
T_K = 298.0
P_PA = 101325.0


def select_phases_al_mg_li(dbf: Database) -> List[str]:
    """Return phases that are compatible with the AL–MG–LI–VA component set.

    We filter to phases whose constituent set is a subset of the allowed
    components in order to keep the equilibrium calculation efficient.
    """
    allowed = {"AL", "MG", "LI", "VA"}
    sel = []
    for ph, pobj in dbf.phases.items():
        const = pobj.constituents
        # Flatten the constituent symbols across sublattices
        flat = set()
        for subl in const:
            for c in subl:
                # Constituents may be Species objects or strings; normalize to symbol/name
                sym = getattr(c, "name", getattr(c, "symbol", str(c))).upper()
                flat.add(sym)
        if flat.issubset(allowed):
            sel.append(ph)
    # A minimal sanity check to ensure we didn't filter out everything
    if not sel:
        raise RuntimeError("No phases compatible with {AL,MG,LI,VA} were found in the database.")
    return sorted(sel)


def stable_phases_summary(eq_ds) -> List[str]:
    """List stable phases at the current equilibrium state (threshold 1e-3).

    The equilibrium Dataset may have a 'vertex' dimension when multiple tie-line
    vertices are present; we sum the per-vertex phase amounts 'NP' to get phase
    moles, then normalize by the total to get phase fractions.
    """
    if "NP" not in eq_ds:
        return []
    pf = eq_ds["NP"]
    if "vertex" in pf.dims:
        pf = pf.sum("vertex")
    # If there is no 'phase' dimension (single phase without labeling), return empty
    if "phase" not in pf.dims:
        return []
    # Sum over phase dimension to get total moles
    tot = pf.sum("phase")
    # Avoid divide-by-zero; if tot is zero, return empty set
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = pf / tot
    # Reduce over all dims except 'phase'
    # We select the first index along any remaining dims as we compute pointwise
    other_dims = [d for d in frac.dims if d != "phase"]
    if other_dims:
        frac0 = frac.isel(**{d: 0 for d in other_dims})
    else:
        frac0 = frac
    stable = (frac0 > 1e-3)
    phases = list(frac0["phase"].where(stable, drop=True).values)
    # Decode bytes to str if needed
    phases = [p.decode() if hasattr(p, "decode") else str(p) for p in phases]
    return phases


def compute_mu_li_alloy(
    dbf: Database, phases: List[str], x_li: float
) -> Dict[str, Any]:
    """Compute mu_Li in Al–Mg–Li at given global X(LI)=x_li with Al:Mg=1:1.

    Returns a dict with keys: x_li, x_al, mu_li_alloy (J/mol), phases (stable list).
    """
    # Enforce composition constraints: for 3 substitutional components, set 2 X()'s
    x_al = 0.5 * (1.0 - x_li)
    conds = {
        v.T: T_K,
        v.P: P_PA,
        v.N: 1.0,
        v.X("LI"): float(x_li),
        v.X("AL"): float(x_al),
    }
    # DOF satisfied by setting X(LI) and X(AL); X(MG) implied.

    comps = ["AL", "MG", "LI", "VA"]
    eq = equilibrium(dbf, comps, phases, conds, verbose=False)
    mu_li = eq.MU.sel(component="LI").values
    # Average across any residual tie-line/vertex indices; this should be uniform
    mu_li_val = float(np.nanmean(mu_li))
    phs = stable_phases_summary(eq)
    return {"x_li": float(x_li), "x_al": float(x_al), "mu_li_alloy": mu_li_val, "phases": phs}


def compute_mu_li_metal(dbf: Database, phases: List[str]) -> float:
    """Compute mu_Li in pure Li metal at 298 K.

    Use the same component set, but constrain X(LI)=1, X(AL)=0 (MG implied 0).
    """
    # Use a reduced component set strictly for Li metal
    comps = ["LI", "VA"]
    # For pure Li, do not set composition constraints (N−1 = 0 for a single substitutional component)
    conds = {v.T: T_K, v.P: P_PA, v.N: 1.0}
    # Restrict to canonical Li-metal phases
    li_phases = [p for p in dbf.phases.keys() if p in {"BCC_A2", "LIQUID"}]
    if not li_phases:
        raise RuntimeError("No Li-metal phases (BCC_A2/LIQUID) found in database.")
    eq = equilibrium(dbf, comps, li_phases, conds, verbose=False)
    mu_li = eq.MU.sel(component="LI").values
    mu_li_val = float(np.nanmean(mu_li))
    if not np.isfinite(mu_li_val):
        raise RuntimeError("Unable to compute mu_Li in Li metal at 298 K (nan).")
    return mu_li_val


def main() -> None:
    # Locate database (cached path; do not re-read CALPHAD.md repeatedly)
    tdb_path = Path("/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB").resolve()
    if not tdb_path.exists():
        raise FileNotFoundError(f"TDB not found: {tdb_path}")

    # Load database once
    dbf = Database(str(tdb_path))

    # Select phases compatible with AL–MG–LI–VA
    phases = select_phases_al_mg_li(dbf)

    # Compute mu_Li in Li metal reference (once)
    mu_li_metal = compute_mu_li_metal(dbf, phases)

    # Composition grid: small x→0 (dilute) to modest Li contents (keep host ratio)
    x_grid = [1e-6, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20]

    rows = []
    for x in x_grid:
        rec = compute_mu_li_alloy(dbf, phases, x)
        # OCV vs Li/Li+: E = - (mu_alloy - mu_metal) / F
        E = - (rec["mu_li_alloy"] - mu_li_metal) / FARADAY
        rec.update({"mu_li_metal": mu_li_metal, "E_V": float(E)})
        rows.append(rec)

    # Prepare summary outputs
    dilute = rows[0]
    summary = {
        "tdb": str(tdb_path),
        "T_K": T_K,
        "P_Pa": P_PA,
        "host_ratio": "Al:Mg = 1:1",
        "definition": "E = - (mu_Li(alloy) - mu_Li(metal)) / F",
        "F_C_per_mol": FARADAY,
        "mu_Li_metal_J_per_mol": mu_li_metal,
        "dilute_limit": {
            "x_li": dilute["x_li"],
            "x_al": dilute["x_al"],
            "mu_li_alloy_J_per_mol": dilute["mu_li_alloy"],
            "E_V": dilute["E_V"],
            "stable_phases": dilute["phases"],
        },
        "grid": rows,
    }

    # Write results
    out_json = Path("results_al_mg_li_ocv_298K.json")
    out_csv = Path("summary_al_mg_li_ocv_298K.csv")
    out_txt = Path("summary_al_mg_li_ocv_298K.txt")

    with out_json.open("w") as f:
        json.dump(summary, f, indent=2)

    # CSV header and rows
    with out_csv.open("w") as f:
        f.write("x_Li,x_Al,mu_Li_alloy_J_per_mol,mu_Li_metal_J_per_mol,E_V,stable_phases\n")
        for r in rows:
            phases_str = "+".join(r["phases"]) if r.get("phases") else ""
            f.write(
                f"{r['x_li']:.6f},{r['x_al']:.6f},{r['mu_li_alloy']:.6f},{mu_li_metal:.6f},{r['E_V']:.6f},{phases_str}\n"
            )

    # Human-readable concise summary
    with out_txt.open("w") as f:
        f.write(f"Using TDB: {tdb_path}\n")
        f.write(f"T = {T_K} K, P = {P_PA} Pa\n")
        f.write("Host ratio fixed: Al:Mg = 1:1 (X(AL)=X(MG) before insertion)\n")
        f.write(f"mu_Li (metal) = {mu_li_metal:.3f} J/mol\n")
        f.write("\nDilute-limit (x_Li → 0) OCV estimate:\n")
        f.write(
            f"  x_Li={dilute['x_li']}, E = {dilute['E_V']:.4f} V, phases: {', '.join(dilute['phases'])}\n"
        )
        f.write("\nGrid results (x_Li, E [V], dominant phases):\n")
        for r in rows:
            f.write(
                f"  x={r['x_li']:.3f}, E={r['E_V']:.4f} V, phases={'+'.join(r['phases']) if r['phases'] else ''}\n"
            )

    # Print concise runtime summary
    print(
        f"Dilute-limit OCV at 298 K (Al:Mg=1:1): {dilute['E_V']:.4f} V; "
        f"mu_Li(metal)={mu_li_metal:.1f} J/mol"
    )
    print(f"Wrote: {out_json}, {out_csv}, {out_txt}")


if __name__ == "__main__":
    main()
