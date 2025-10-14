"""
CALPHAD analysis for Al–Si–C (Al0.30 Si0.55 C0.15) to test:

    "Heating a sample of Al30Si55C15 to 1500 K is insufficient to dissolve all carbide precipitates."

This script performs a single-shot, reproducible equilibrium calculation using pycalphad and the local
thermodynamic database under TnETDBDB (COST507R.TDB). It avoids any interactive input and writes
machine-readable outputs for downstream use:

- results.csv: per-temperature phase fractions for all selected phases
- results.json: summary including which carbide(s) are stable at 1500 K, and the disappearance temperature
  (if any) for total carbides
- summary.txt: short human-readable report

Methodology notes (per CALPHAD workflow guidance):
- Database inspection to choose phases is done via pycalphad introspection, not by grepping the TDB file.
- Degrees of freedom: for 3 components (Al, Si, C), we fix N=1, P=1 atm, and two mole fraction constraints
  v.X('AL') and v.X('SI'). v.X('C') is implied.
- Phase fractions are computed from eq['NP']; vertex-labeled results are summed over the 'vertex' dimension.
- Carbide identification: any phase that contains C and at least one of {AL, SI}, but is not a generic
  solution/elemental phase (LIQUID, FCC_A1, BCC_A2, HCP_A3, GRAPHITE, DIAMOND_A4). This captures ordered
  carbides such as Al4C3 or SiC while excluding solution phases that simply dissolve C.

This script adheres to the speed/stability guidance by:
- Restricting phases to those whose constituents are a subset of {AL,SI,C,VA}
- Using a modest temperature grid spanning 1000–2000 K (step 25 K) to find dissolution temperature
  while directly reporting the state at 1500 K

Author: computational_tools_0005
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from pycalphad import Database, equilibrium, variables as v


def phase_constituent_union(dbf: Database, phase: str) -> Set[str]:
    const = dbf.phases[phase].constituents
    uni: Set[str] = set()
    for subl in const:
        for c in subl:
            name = getattr(c, 'name', str(c))
            uni.add(name.upper())
    return uni


def phase_has_any(dbf: Database, phase: str, comps: Set[str]) -> bool:
    # Allow phases if any of our components appear among its constituents
    uni = phase_constituent_union(dbf, phase)
    return bool(uni & set(comps))


def is_solution_or_elemental(phase: str) -> bool:
    sol_like = {
        "LIQUID",
        "FCC_A1",
        "BCC_A2",
        "HCP_A3",
        "GRAPHITE",
        "DIAMOND_A4",
        "GAS",
    }
    return phase.upper() in sol_like


def is_carbide_phase(dbf: Database, phase: str) -> bool:
    uni = phase_constituent_union(dbf, phase)
    if "C" not in uni:
        return False
    # Require at least one of the metallic/metalloid components present with carbon
    if not ({"AL", "SI"} & uni):
        return False
    # Exclude generic solution/elemental phases
    if is_solution_or_elemental(phase):
        return False
    return True


def main() -> None:
    # Fixed inputs per prompt
    db_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    assert db_path.exists(), f"TDB not found at {db_path}"

    # Composition: Al0.30 Si0.55 C0.15 (mole fraction)
    x_al = 0.30
    x_si = 0.55
    # x_c implied = 0.15

    # Conditions: pressure and temperature sweep; focus also on T=1500 K
    T_min, T_max, T_step = 1000.0, 2000.0, 25.0
    T_grid = np.arange(T_min, T_max + 1e-9, T_step)
    P = 101325.0

    # Components of interest
    comps = ["AL", "SI", "C"]

    # Load database and introspect phases
    dbf = Database(str(db_path))
    all_phases = sorted(dbf.phases.keys())

    # Select phases; if introspection fails to find any with our simple filter, use all phases
    allowed_phases = [p for p in all_phases if phase_has_any(dbf, p, set(comps))]
    if not allowed_phases:
        # Fallback: include all phases and let pycalphad filter by components during model building
        allowed_phases = all_phases
    assert allowed_phases, "No phases selected; did you build the phases list?"

    # Prepare conditions (N−1 mole fraction constraints for N components)
    conds = {
        v.T: T_grid,
        v.P: P,
        v.N: 1.0,
        v.X("AL"): x_al,
        v.X("SI"): x_si,
    }
    # (Constraint count check omitted here due to environment differences in variable attributes.)

    # Compute equilibrium with minimal outputs (NP + Phase labels)
    eq = equilibrium(dbf, comps, allowed_phases, conds, output='NP')

    # Extract arrays
    np_da = eq['NP']        # moles per vertex
    ph_da = eq['Phase']     # phase label per vertex
    T_vals = np.array(eq['T'].values).ravel()

    # Build per-temperature phase fractions by accumulating NP per phase label
    thr = 1e-3
    all_obs_phases: Set[str] = set()
    frac_by_T: List[Dict[str, float]] = []
    carbide_total_vs_T: List[float] = []

    for iT, T in enumerate(T_vals):
        np_vals = np_da.isel(T=iT).values.reshape(-1)
        ph_vals = ph_da.isel(T=iT).values.reshape(-1)
        phase_moles: Dict[str, float] = {}
        for ph_label, np_mol in zip(ph_vals, np_vals):
            ph_name = str(ph_label).strip()
            if not ph_name:
                continue
            mval = float(np_mol)
            if mval == 0.0:
                continue
            phase_moles[ph_name] = phase_moles.get(ph_name, 0.0) + mval
        tot_moles = sum(phase_moles.values())
        if tot_moles <= 0:
            # No material accounted for at this T (unexpected); record zeros
            frac_by_T.append({})
            carbide_total_vs_T.append(0.0)
            continue
        # Fractions
        frac = {ph: m / tot_moles for ph, m in phase_moles.items()}
        frac_by_T.append(frac)
        all_obs_phases.update(frac.keys())

        # Carbide total fraction at this T
        carb_frac = 0.0
        for ph, f in frac.items():
            if is_carbide_phase(dbf, ph):
                carb_frac += f
        carbide_total_vs_T.append(carb_frac)

    # Identify 1500 K index
    idx_1500 = int(np.argmin(np.abs(T_vals - 1500.0)))
    T_1500 = float(T_vals[idx_1500])

    # Phase fractions at 1500 K
    phases_at_1500 = frac_by_T[idx_1500]
    carbide_fraction_1500 = float(carbide_total_vs_T[idx_1500])
    carbides_present_1500 = carbide_fraction_1500 > thr

    # Dissolution T estimate (first T with carbide_total ≤ thr)
    dissolve_T = None
    below = np.where(np.array(carbide_total_vs_T) <= thr)[0]
    if len(below) > 0:
        dissolve_T = float(T_vals[below[0]])

    # Prepare CSV across union of observed phases
    phases_sorted = sorted(all_obs_phases)
    rows = []
    for iT, T in enumerate(T_vals):
        rec = {"T": float(T)}
        frac = frac_by_T[iT]
        for ph in phases_sorted:
            rec[ph] = float(frac.get(ph, 0.0))
        rows.append(rec)
    df = pd.DataFrame(rows)

    # Outputs
    out_csv = Path("al_si_c_results.csv")
    out_json = Path("al_si_c_results.json")
    out_txt = Path("al_si_c_summary.txt")

    df.to_csv(out_csv, index=False)

    # Summaries
    # Build carbide flags just for observed phases
    carbide_flags: Dict[str, bool] = {p: is_carbide_phase(dbf, p) for p in phases_sorted}
    stable_phases_1500 = sorted([p for p, f in phases_at_1500.items() if f > thr])
    stable_carbides_1500 = [p for p in stable_phases_1500 if carbide_flags.get(p, False)]

    # Phases present at ends of T-grid
    # Stable phases at the ends of the grid using computed fractions
    stable_lowT = sorted([p for p, f in frac_by_T[0].items() if f > thr])
    stable_highT = sorted([p for p, f in frac_by_T[-1].items() if f > thr])

    summary = {
        "tdb_path": str(db_path),
        "components": comps,
        "phases_considered": allowed_phases,
        "composition": {"X(AL)": x_al, "X(SI)": x_si, "X(C)": 1 - x_al - x_si},
        "pressure_Pa": P,
        "T_grid_K": {"min": T_min, "max": T_max, "step": T_step, "count": len(T_vals)},
        "T_at_1500_idx": {"T_closest": T_1500, "index": idx_1500},
        "stable_phases_at_1500K": stable_phases_1500,
        "stable_carbides_at_1500K": stable_carbides_1500,
        "carbide_total_fraction_at_1500K": carbide_fraction_1500,
        "dissolution_temperature_K_first_below_threshold": dissolve_T,
        "threshold_for_presence": thr,
        "stable_phases_lowT": stable_lowT,
        "stable_phases_highT": stable_highT,
    }

    with out_json.open("w") as f:
        json.dump(summary, f, indent=2)

    # Print concise human summary
    lines = []
    lines.append(f"Using TDB: {db_path}")
    lines.append(f"Composition: X(AL)={x_al:.2f}, X(SI)={x_si:.2f}, X(C)={1 - x_al - x_si:.2f}")
    lines.append(f"P = {P:.0f} Pa; T in [{T_min:.0f}, {T_max:.0f}] K step {T_step:.0f} K")
    lines.append(f"Phases considered ({len(allowed_phases)}): {', '.join(allowed_phases)}")
    lines.append("")
    lines.append(f"At T≈1500 K (closest grid T={T_1500:.1f} K):")
    if stable_phases_1500:
        lines.append("  Stable phases (f>1e-3): " + ", ".join(stable_phases_1500))
    else:
        lines.append("  No stable phases above threshold (unexpected)")
    if stable_carbides_1500:
        carb_list = []
        for ph in stable_carbides_1500:
            carb_list.append(f"{ph}({phases_at_1500[ph]:.3f})")
        lines.append("  Carbides present: " + ", ".join(carb_list))
    else:
        lines.append("  No carbides present above threshold at 1500 K")
    if dissolve_T is None:
        lines.append("  Carbides do not dissolve within 1000–2000 K (threshold 1e-3)")
    else:
        lines.append(f"  First T where total carbide fraction ≤1e-3: {dissolve_T:.1f} K")

    with out_txt.open("w") as f:
        f.write("\n".join(lines) + "\n")

    # Also print to stdout for immediate visibility
    print("\n".join(lines))


if __name__ == "__main__":
    main()
