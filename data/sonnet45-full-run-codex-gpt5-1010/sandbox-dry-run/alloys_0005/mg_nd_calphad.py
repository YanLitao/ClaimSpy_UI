#!/usr/bin/env python
"""
Mg–Nd CALPHAD analysis using a local TDB (COST507R) to determine:
  1) Stable phases and phase fractions at 40 wt.% Nd, 298 K
  2) Density at that composition and temperature (approximate via ideal volume mixing)
  3) Stable phases and fractions at 40 at.% Nd, 298 K (comparison)

Outputs
  - results.json
  - results.csv
  - Prints a concise human summary to stdout

Notes on density: The COST507R database does not expose an immediately
accessible molar-volume property across all phases via pycalphad. Here we
estimate alloy density by ideal volume mixing of the pure elements at 298 K,
independent of phase constitution. This yields a conservative, reproducible
estimate:
    ρ_mix = (total mass) / (Σ_i mass_i / ρ_i)
where i ∈ {Mg, Nd} with masses determined from the requested composition
definition (wt.% or at.%).

This script loads the TDB exactly once and solves both requested conditions
in one pass, following the CALPHAD.md workflow.
"""

from __future__ import annotations

import json
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from pycalphad import Database, equilibrium, variables as v


# -------------------------
# Configuration and helpers
# -------------------------

# Absolute path to local TDB as required by CALPHAD.md
TDB_PATH = "/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB"

# Components and pressure
COMPS = ["MG", "ND", "VA"]
PRESSURE = 101_325  # Pa

# Threshold for reporting a phase as present
PHASE_PRESENCE_TOL = 1e-3

# Pure-element data at ~298 K (room temperature)
# Densities from standard references (g/cm^3); atomic weights (g/mol)
AW = {"MG": 24.305, "ND": 144.242}
RHO_298K = {"MG": 1.738, "ND": 7.007}


def wt_to_molefrac(wt_nd: float, aw_mg=AW["MG"], aw_nd=AW["ND"]) -> float:
    """Convert Nd weight fraction to Nd mole fraction in Mg–Nd binary.
    wt_nd is in [0,1].
    """
    wt_mg = 1.0 - wt_nd
    n_nd = wt_nd / aw_nd
    n_mg = wt_mg / aw_mg
    return n_nd / (n_nd + n_mg)


def ideal_mixture_density_from_wt(wt_nd: float) -> float:
    """Density (g/cm^3) by ideal volume mixing for Mg–Nd at 298 K given Nd weight fraction.
    """
    wt_mg = 1.0 - wt_nd
    mass_total = 100.0  # g, arbitrary scale
    m_nd = mass_total * wt_nd
    m_mg = mass_total * wt_mg
    vol_total_cm3 = m_nd / RHO_298K["ND"] + m_mg / RHO_298K["MG"]
    return mass_total / vol_total_cm3


def ideal_mixture_density_from_at(x_nd: float) -> float:
    """Density (g/cm^3) by ideal volume mixing for Mg–Nd at 298 K given Nd atomic fraction.
    """
    x_mg = 1.0 - x_nd
    # Work per 1 mole of atoms
    m_nd = x_nd * AW["ND"]  # g
    m_mg = x_mg * AW["MG"]  # g
    vol_total_cm3 = m_nd / RHO_298K["ND"] + m_mg / RHO_298K["MG"]
    mass_total = m_nd + m_mg
    return mass_total / vol_total_cm3


def build_phase_list(db: Database, comps: List[str]) -> List[str]:
    """Select phases from the database that contain at least one of the components.
    Following CALPHAD.md guidance to avoid grepping files—use pycalphad introspection.
    """
    def phase_has_any(phase: str, comps_: List[str]) -> bool:
        const = db.phases[phase].constituents
        # Normalize to string symbols for comparison
        flat = {str(c).strip().upper() for subl in const for c in subl}
        return bool(flat & {c.upper() for c in comps_})

    phases = sorted(p for p in db.phases.keys() if phase_has_any(p, ["MG", "ND"]))
    if not phases:
        raise RuntimeError("No phases selected; did you build the phases list?")
    return phases


def summarize_equilibrium(eq, temp_K: float) -> Dict[str, float]:
    """Return phase fractions as a dict {phase: fraction} for a single (T,P,x) point.
    Handles modern pycalphad outputs where phases are enumerated over the 'vertex' axis.
    """
    # Extract along single-condition indices
    pf = eq["NP"]
    # Reduce singleton dims N,P,T,X_* to scalars, keep 'vertex'
    for dim in [d for d in pf.dims if d != "vertex"]:
        if pf.sizes.get(dim, 1) == 1:
            pf = pf.isel({dim: 0})

    # Corresponding phases per vertex
    phases_da = eq["Phase"]
    for dim in [d for d in phases_da.dims if d != "vertex"]:
        if phases_da.sizes.get(dim, 1) == 1:
            phases_da = phases_da.isel({dim: 0})

    # Collect NP per phase over vertices
    def _norm_name(p):
        return p.decode() if isinstance(p, (bytes, bytearray)) else str(p)

    per_phase_moles: Dict[str, float] = {}
    total_moles = 0.0
    for iv in range(pf.sizes.get("vertex", 0)):
        name = _norm_name(phases_da.values[iv])
        if not name:
            continue
        n = float(pf.values[iv])
        if not np.isfinite(n):
            continue
        total_moles += n
        per_phase_moles[name] = per_phase_moles.get(name, 0.0) + n

    if total_moles <= 0:
        return {}

    # Normalize to fractions and threshold
    present = {
        ph: n / total_moles for ph, n in per_phase_moles.items() if (n / total_moles) > PHASE_PRESENCE_TOL
    }
    return dict(sorted(present.items(), key=lambda kv: kv[1], reverse=True))


def main():
    out_dir = Path(".")
    tdb_path = Path(TDB_PATH).resolve()
    assert tdb_path.exists(), f"TDB not found: {tdb_path}"

    # Load DB once
    db = Database(str(tdb_path))
    phases = build_phase_list(db, COMPS)

    # Condition set for equilibrium
    base_conds = {v.P: PRESSURE, v.N: 1.0}

    # Requested conditions
    T_room = 298.0

    # Case A: 40 wt.% Nd → mole fraction
    wt_nd = 0.40
    x_nd_from_wt = wt_to_molefrac(wt_nd)

    conds_wt = dict(base_conds)
    conds_wt.update({v.T: T_room, v.X("ND"): x_nd_from_wt})

    # Case B: 40 at.% Nd
    x_nd = 0.40
    conds_at = dict(base_conds)
    conds_at.update({v.T: T_room, v.X("ND"): x_nd})

    # Solve both equilibria (separately) to keep memory bounded
    eq_wt = equilibrium(db, COMPS, phases, conds_wt, calc_opts={"pdens": 1200})
    eq_at = equilibrium(db, COMPS, phases, conds_at, calc_opts={"pdens": 1200})

    # Summaries
    # Summaries
    try:
        phases_wt = summarize_equilibrium(eq_wt, T_room)
    except Exception as e:
        # Fallback: derive presence from mask if reduction across 'phase' fails
        mask = (eq_wt["NP"] > PHASE_PRESENCE_TOL)
        phvals = eq_wt["Phase"].where(mask, drop=True).values
        uniq = sorted(set([p.decode() if isinstance(p, (bytes, bytearray)) else str(p) for p in np.ravel(phvals)]))
        phases_wt = {p: float('nan') for p in uniq}

    try:
        phases_at = summarize_equilibrium(eq_at, T_room)
    except Exception as e:
        mask = (eq_at["NP"] > PHASE_PRESENCE_TOL)
        phvals = eq_at["Phase"].where(mask, drop=True).values
        uniq = sorted(set([p.decode() if isinstance(p, (bytes, bytearray)) else str(p) for p in np.ravel(phvals)]))
        phases_at = {p: float('nan') for p in uniq}

    # Densities via ideal volume mixing at 298 K
    rho_wt = ideal_mixture_density_from_wt(wt_nd)
    rho_at = ideal_mixture_density_from_at(x_nd)

    # Write CSV and JSON
    rows = []
    for label, xnd, phases_map, rho in [
        ("40 wt% Nd", x_nd_from_wt, phases_wt, rho_wt),
        ("40 at% Nd", x_nd, phases_at, rho_at),
    ]:
        # One row per phase for clarity
        if phases_map:
            for ph, frac in phases_map.items():
                rows.append({
                    "case": label,
                    "T_K": T_room,
                    "P_Pa": PRESSURE,
                    "x_Nd_mole": xnd,
                    "phase": ph,
                    "phase_fraction": frac,
                    "density_g_cm3": rho,
                    "tdb_path": str(tdb_path),
                })
        else:
            rows.append({
                "case": label,
                "T_K": T_room,
                "P_Pa": PRESSURE,
                "x_Nd_mole": xnd,
                "phase": "(none)",
                "phase_fraction": 0.0,
                "density_g_cm3": rho,
                "tdb_path": str(tdb_path),
            })

    csv_path = out_dir / "results.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "case", "T_K", "P_Pa", "x_Nd_mole", "phase", "phase_fraction", "density_g_cm3", "tdb_path"
            ],
        )
        w.writeheader()
        w.writerows(rows)

    json_obj = {
        "tdb": {
            "chosen": str(tdb_path),
            "considered": [str(tdb_path)],
        },
        "conditions": {
            "T_K": T_room,
            "P_Pa": PRESSURE,
        },
        "cases": [
            {
                "label": "40 wt% Nd",
                "x_Nd_mole": x_nd_from_wt,
                "phases": phases_wt,
                "density_g_cm3": rho_wt,
            },
            {
                "label": "40 at% Nd",
                "x_Nd_mole": x_nd,
                "phases": phases_at,
                "density_g_cm3": rho_at,
            },
        ],
        "assumptions": {
            "pressure": "Assumed 101325 Pa (1 atm)",
            "density_model": "Ideal volume mixing of pure elements at 298 K",
            "phase_threshold": PHASE_PRESENCE_TOL,
        },
    }
    json_path = out_dir / "results.json"
    with json_path.open("w") as f:
        json.dump(json_obj, f, indent=2)

    # Human summary to stdout
    def fmt_phases(d: Dict[str, float]) -> str:
        if not d:
            return "(no stable phases above threshold)"
        return ", ".join([f"{k}: {v:.3f}" for k, v in d.items()])

    print("Database:", tdb_path)
    print(f"Conditions: T={T_room:.1f} K, P={PRESSURE} Pa")
    print(f"40 wt.% Nd → x_Nd = {x_nd_from_wt:.5f}")
    print("  Phases:", fmt_phases(phases_wt))
    print(f"  Density (ideal mix, 298 K): {rho_wt:.4f} g/cm^3")
    print(f"40 at.% Nd → x_Nd = {x_nd:.5f}")
    print("  Phases:", fmt_phases(phases_at))
    print(f"  Density (ideal mix, 298 K): {rho_at:.4f} g/cm^3")
    print("Wrote:", json_path)
    print("Wrote:", csv_path)


if __name__ == "__main__":
    main()
