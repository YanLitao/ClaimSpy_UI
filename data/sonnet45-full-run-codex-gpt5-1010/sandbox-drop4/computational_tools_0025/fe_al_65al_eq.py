import json
from pathlib import Path

import numpy as np
import pandas as pd
from pycalphad import Database, equilibrium
from pycalphad import variables as v


def main():
    # Cache essentials from CALPHAD.md: database dir and env
    tdb_path = "/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB"

    # System definition
    components = ["FE", "AL", "VA"]
    x_al = 0.65  # 65 at.% Al
    temps = [298.0, 900.0]

    # Load database once
    db = Database(tdb_path)

    # Build a relevant phase list: phases that have any of FE/AL in their constituents
    # (solution phases often list many elements; intersection is appropriate)
    target = {"FE", "AL"}

    def phase_has_any(phase_name):
        const = db.phases[phase_name].constituents
        present = {getattr(el, 'name', str(el)) for subl in const for el in subl}
        return bool(present & target)

    phases = sorted(p for p in db.phases if phase_has_any(p))

    # Fail fast if nothing selected
    if not phases:
        print(f"Loaded database from: {tdb_path}")
        print(f"Total phases in DB: {len(db.phases)}")
        # Show a few phase names to aid debugging
        try:
            print("Sample phases:", list(db.phases.keys())[:10])
        except Exception:
            pass
        raise AssertionError("No phases selected; verify the Fe–Al database path and contents.")

    # Equilibrium conditions (fix N−1 site fractions for a binary)
    conds = {v.T: temps, v.P: 101325.0, v.N: 1.0, v.X("AL"): x_al}

    # Compute equilibrium once for both temperatures
    eq = equilibrium(db, components, phases, conds)
    # Debug: show NP dims if phase dimension is missing
    try:
        np_da = eq["NP"]
        print("NP dims:", np_da.dims)
        print("Coordinates:", list(eq.coords))
        print("Data variables:", list(eq.data_vars))
        if "Phase" in eq:
            print("Phase dims:", eq["Phase"].dims)
    except Exception:
        pass

    # Helper to squeeze all non-vertex dims to their first index
    def squeeze_to_vertex(da):
        out = da
        for d in list(out.dims):
            if d != "vertex":
                out = out.isel({d: 0})
        return out

    # Collect per‑temperature results
    results = []
    summary_by_T = {}
    for iT, T in enumerate(temps):
        # Build per-phase totals by accumulating NP over vertices
        phase_var = eq["Phase"].sel(T=T)
        np_var = eq["NP"].sel(T=T)
        phase_1d = squeeze_to_vertex(phase_var).values
        np_1d = squeeze_to_vertex(np_var).values

        # Ensure 1D arrays over 'vertex'
        phase_1d = np.ravel(phase_1d)
        np_1d = np.ravel(np_1d)

        totals = {}
        for name, amount in zip(phase_1d, np_1d):
            if not np.isfinite(amount) or amount <= 0:
                continue
            # phase names may be bytes/Species; convert to upper-case string
            pname = str(getattr(name, 'name', name)).upper()
            totals[pname] = totals.get(pname, 0.0) + float(amount)

        tot_moles = sum(totals.values())
        stable = []
        if tot_moles > 0:
            for p, n in totals.items():
                frac = n / tot_moles
                if frac > 1e-3:
                    stable.append((p, frac))
        stable.sort(key=lambda x: -x[1])

        # Identify BCC-based primary (A2 or B2)
        primary_phase = stable[0][0] if stable else None
        bcc_primary = primary_phase in {"BCC_A2", "B2_BCC", "BCC_B2"}

        # Secondary Al2Fe presence
        al2fe_present = any(p.startswith("AL2FE") or p == "AL2FE" for p, _ in stable)

        summary_by_T[T] = {
            "primary_phase": primary_phase,
            "is_bcc_primary": bool(bcc_primary),
            "al2fe_present": bool(al2fe_present),
            "phases": [{"phase": p, "mole_fraction": f} for p, f in stable],
        }

        for p, f in stable:
            results.append({"T": T, "phase": p, "mole_fraction": f})

    # Write outputs
    out_json = Path("fe_al_65al_results.json")
    out_csv = Path("fe_al_65al_results.csv")
    with out_json.open("w") as fp:
        json.dump({
            "database": str(tdb_path),
            "composition": {"X_AL": x_al, "X_FE": 1 - x_al},
            "pressure_Pa": 101325.0,
            "temperatures_K": temps,
            "results_by_temperature": summary_by_T,
        }, fp, indent=2)
    pd.DataFrame(results).to_csv(out_csv, index=False)

    # Print concise human summary
    print("Fe–Al (X_Al = 0.65) equilibrium summary")
    print(f"Database: {tdb_path}")
    for T in temps:
        s = summary_by_T[T]
        print(f"\nT = {T:.0f} K @ 1 atm")
        print(f"  Primary phase: {s['primary_phase']}")
        print(f"  BCC-based primary (A2/B2): {s['is_bcc_primary']}")
        print(f"  Al2Fe present: {s['al2fe_present']}")
        for item in s["phases"]:
            print(f"    - {item['phase']}: mole fraction = {item['mole_fraction']:.4f}")


if __name__ == "__main__":
    main()
