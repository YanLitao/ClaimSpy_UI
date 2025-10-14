import os
import json
import numpy as np
import pandas as pd

from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf, phase, comps):
    const = dbf.phases[phase].constituents
    flat = {c for subl in const for c in subl}
    return bool(set(comps) & flat)


def compute_phase_fractions(eq):
    pf = eq["NP"]
    if "vertex" in pf.dims:
        pf = pf.sum("vertex")
    tot = pf.sum("phase")
    # Avoid division by zero
    phase_frac = pf / tot.where(tot != 0)
    return phase_frac


def structure_label(phase_name: str) -> str:
    name = phase_name.upper()
    if "LIQUID" in name:
        return "LIQUID"
    if "FCC" in name:
        return "FCC"
    if "HCP" in name:
        return "HCP"
    if "BCC" in name or "A2" in name:
        return "BCC"
    if "L12" in name:
        return "FCC-ordered (L12)"
    if "L10" in name:
        return "FCC-ordered (L10)"
    if "OMEGA" in name:
        return "HCP-derivative (omega)"
    return "Intermetallic/other"


def summarize_at_T(eq, T_value, threshold=1e-3):
    # Select nearest T index
    T_coords = eq.coords["T"].values
    idx = int(np.argmin(np.abs(T_coords - T_value)))
    T_actual = float(T_coords[idx])
    phase_frac = compute_phase_fractions(eq)
    # Reduce across all non-phase dims except 'phase' and 'T'
    reduce_dims = [d for d in phase_frac.dims if d not in ("phase", "T")]
    pf_T = phase_frac.isel(T=idx)
    if reduce_dims:
        pf_T = pf_T.mean(dim=reduce_dims)
    # Build summary dictionary of phase fractions
    result = {}
    for ph in pf_T.coords["phase"].values:
        val = float(pf_T.sel(phase=ph).values)
        if np.isfinite(val) and val > threshold:
            result[str(ph)] = val
    # Sort by fraction descending
    result = dict(sorted(result.items(), key=lambda kv: kv[1], reverse=True))
    return T_actual, result


def main():
    # Inputs
    workdir = os.getcwd()
    out_json = os.path.join(workdir, "results.json")
    out_csv = os.path.join(workdir, "results.csv")
    out_txt = os.path.join(workdir, "summary.txt")

    # Composition (atomic fraction)
    comps = ["AL", "LI", "MG", "SC", "TI"]
    x_spec = {"AL": 0.20, "LI": 0.20, "MG": 0.10, "SC": 0.20}  # TI implied (0.30)

    # Database path (must be local TnETDBDB)
    tdb_path = "/Users/delip/play/k2code/TnETDBDB/cost507R.TDB"

    summary_lines = []
    results = {
        "tdb": tdb_path,
        "composition_atpct": {"AL": 20, "LI": 20, "MG": 10, "SC": 20, "TI": 30},
        "pressure_Pa": 101325,
        "notes": [],
        "conditions": {},
    }

    if not os.path.exists(tdb_path):
        msg = f"TDB not found: {tdb_path}"
        results["error"] = msg
        with open(out_json, "w") as f:
            json.dump(results, f, indent=2)
        print(msg)
        return

    # Load DB once
    dbf = Database(tdb_path)
    db_elems = sorted({str(e).upper() for e in dbf.elements})
    missing = [el for el in comps if el not in db_elems]
    results["db_elements"] = db_elems
    if missing:
        msg = (
            "Selected TDB does not cover required elements; cannot evaluate this alloy. "
            f"Missing: {missing}"
        )
        results["error"] = msg
        results["notes"].append(
            "Per workflow: Only TDBs under TnETDBDB may be used; this one lacks required components."
        )
        with open(out_json, "w") as f:
            json.dump(results, f, indent=2)
        print(msg)
        return

    # Build phases list filtered by constituents to keep the problem tractable
    all_phases = sorted(dbf.phases.keys())
    phases = [p for p in all_phases if phase_has_any(dbf, p, comps)]
    assert phases, "No phases selected; did you build the phases list?"

    # Conditions: Fix N−1 mole fractions
    conds_common = {v.P: 101325, v.N: 1.0}
    conds_common.update({v.X(k): v for k, v in x_spec.items()})

    # Temperature grid includes the requested anneal temperature explicitly
    T_list = np.unique(np.array(list(np.linspace(2000, 500, 31)) + [773.0], dtype=float))

    # Run a single equilibrium sweep across T to capture liquidus vicinity and the 773 K state
    eq = equilibrium(dbf, comps, phases, {**conds_common, v.T: T_list}, verbose=False)

    # High-temperature near-liquidus assessment
    phase_frac = compute_phase_fractions(eq)
    # Identify the highest T where any solid phase fraction is present (>1e-3)
    T_coords = eq.coords["T"].values
    highT_idx = None
    highT_phases = None
    for i in range(len(T_coords) - 1, -1, -1):
        Tval = float(T_coords[i])
        T_actual, pf_dict = summarize_at_T(eq, Tval)
        # Count non-liquid phases
        solids = {k: v for k, v in pf_dict.items() if "LIQUID" not in k.upper()}
        if solids:
            highT_idx = i
            highT_phases = {"T": T_actual, "phases": solids}
            break
    cond1 = {}
    if highT_idx is not None and highT_phases is not None:
        # Determine dominant solid just below liquidus
        dominant_phase, dominant_frac = next(iter(highT_phases["phases"].items()))
        cond1 = {
            "temperature_K": highT_phases["T"],
            "dominant_solid_phase": dominant_phase,
            "dominant_solid_structure": structure_label(dominant_phase),
            "phase_fractions": highT_phases["phases"],
        }
    else:
        cond1 = {"error": "No solid phases found in the scanned T range."}
    results["conditions"]["near_liquidus"] = cond1

    # Equilibrium at anneal temperature 773 K
    T_anneal = 773.0
    T_anneal_actual, pf_anneal = summarize_at_T(eq, T_anneal)
    cond3 = {
        "temperature_K": T_anneal_actual,
        "phase_fractions": pf_anneal,
        "structures": {ph: structure_label(ph) for ph in pf_anneal},
    }
    results["conditions"]["anneal_773K_equilibrium"] = cond3

    # Rapid quench to room temperature (metastable retention of high-T dominant solid)
    T_quench = 298.0
    cond2 = {"temperature_K": T_quench}
    if "dominant_solid_phase" in cond1:
        primary = cond1["dominant_solid_phase"]
        # If the primary phase is not in the phases list (should be), we still try
        phases_metastable = [primary]
        try:
            eq_q = equilibrium(dbf, comps, phases_metastable, {**conds_common, v.T: T_quench}, verbose=False)
            # If calculation succeeds, it implies a metastable single-phase solution exists
            _, pf_q = summarize_at_T(eq_q, T_quench)
            cond2.update({
                "metastable_primary": primary,
                "metastable_primary_structure": structure_label(primary),
                "phase_fractions": pf_q,
                "note": "Metastable calculation with only the primary high-T solid phase allowed. Represents as-quenched retention if diffusion is suppressed.",
            })
        except Exception as e:
            cond2.update({
                "metastable_primary": primary,
                "error": f"Metastable calculation failed: {e}",
            })
    else:
        cond2.update({"error": "No identified primary solid near liquidus to retain upon quench."})
    results["conditions"]["quenched_298K_metastable"] = cond2

    # Build concise human summary
    summary_lines.append(f"Database: {tdb_path}")
    summary_lines.append(f"Elements in DB: {', '.join(results['db_elements'])}")
    if "error" in results:
        summary_lines.append(results["error"])
    else:
        # Condition 1
        if "error" in cond1:
            summary_lines.append("Near-liquidus: could not identify a solid phase in scanned range.")
        else:
            p0, f0 = cond1["dominant_solid_phase"], cond1["dominant_solid_structure"]
            summary_lines.append(
                f"Near-liquidus (~{cond1['temperature_K']:.1f} K): dominant solid {p0} ({f0})."
            )
        # Condition 2
        if "error" in cond2:
            summary_lines.append(f"Quenched 298 K: {cond2['error']}")
        else:
            summary_lines.append(
                f"Quenched 298 K (metastable): retain {cond2['metastable_primary']} ({cond2['metastable_primary_structure']})."
            )
        # Condition 3
        if pf_anneal:
            # Identify dominant phase at 773 K
            dom_anneal = next(iter(sorted(pf_anneal.items(), key=lambda kv: kv[1], reverse=True)))
            summary_lines.append(
                f"Anneal 773 K: dominant phase {dom_anneal[0]} ({structure_label(dom_anneal[0])})."
            )

    # Write outputs
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)

    # Flatten results to CSV rows
    rows = []
    def add_rows(label, T, pf_dict):
        if not pf_dict:
            rows.append({"condition": label, "T_K": T, "phase": "(none)", "fraction": 0.0, "structure": ""})
        else:
            for ph, frac in pf_dict.items():
                rows.append({
                    "condition": label,
                    "T_K": T,
                    "phase": ph,
                    "fraction": frac,
                    "structure": structure_label(ph)
                })

    if "error" not in results:
        if "temperature_K" in cond1:
            add_rows("near_liquidus", cond1["temperature_K"], cond1.get("phase_fractions", {}))
        add_rows("anneal_773K_equilibrium", cond3.get("temperature_K", 773.0), cond3.get("phase_fractions", {}))
        if "phase_fractions" in cond2:
            add_rows("quenched_298K_metastable", cond2["temperature_K"], cond2.get("phase_fractions", {}))

    pd.DataFrame(rows).to_csv(out_csv, index=False)

    with open(out_txt, "w") as f:
        f.write("\n".join(summary_lines) + "\n")

    # Print concise summary to stdout
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()

