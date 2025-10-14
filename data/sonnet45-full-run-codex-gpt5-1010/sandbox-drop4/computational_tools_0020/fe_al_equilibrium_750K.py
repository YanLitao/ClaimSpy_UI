import json
import csv
import sys
from pathlib import Path

import numpy as np

# Fail fast: import pycalphad only once here
try:
    from pycalphad import Database, variables as v, equilibrium
except Exception as e:
    print("ERROR: pycalphad import failed. Ensure the readertools env is active.")
    raise


def find_tdb_dir(start: Path) -> Path:
    cur = start.resolve()
    for _ in range(10):  # climb up to 10 levels
        candidate = cur / 'TnETDBDB'
        if candidate.is_dir():
            return candidate
        if cur.parent == cur:
            break
        cur = cur.parent
    raise FileNotFoundError("TnETDBDB directory not found by climbing parents from: " + str(start))


def phase_has_any(dbf: Database, phase: str, comps: list[str]) -> bool:
    const = dbf.phases[phase].constituents
    flat = {c for subl in const for c in subl}
    return bool(flat & set(comps))


def main():
    cwd = Path.cwd()
    out_json = cwd / 'fe_al_750K_results.json'
    out_csv = cwd / 'fe_al_750K_results.csv'

    # Locate database
    tdb_dir = find_tdb_dir(cwd)
    # Prefer COST507R if present
    tdb_path = tdb_dir / 'cost507R.TDB'
    if not tdb_path.is_file():
        raise FileNotFoundError(f"No cost507R.TDB in {tdb_dir}")

    # Load database
    dbf = Database(str(tdb_path))

    # System definition: Fe–Al binary
    comps = ['AL', 'FE']

    # Build phases list containing our components
    phases = sorted(p for p in dbf.phases if phase_has_any(dbf, p, comps))
    if not phases:
        # Debug information to help diagnose database coverage
        print("DEBUG: No phases matched filter; database elements:", sorted(dbf.elements))
        # Show a few phases and their constituents
        sample = list(sorted(dbf.phases))[:10]
        for ph in sample:
            const = dbf.phases[ph].constituents
            print(f"DEBUG: {ph} constituents: {const}")
        # Fallback to using all phases and let pycalphad handle filtering internally
        phases = sorted(dbf.phases)

    # Conditions at 750 K, 1 atm, X(AL)=0.70 (implies X(Fe)=0.30)
    conds = {v.T: 750.0, v.P: 101325.0, v.N: 1.0, v.X('AL'): 0.70}

    # Compute equilibrium
    # Compute equilibrium; request default outputs (includes NP)
    eq = equilibrium(dbf, comps, phases, conds)

    # Phase moles and fractions
    pf = eq['NP']
    if 'phase' not in pf.dims:
        # Newer pycalphad: NP is per-vertex; phase identity is in eq['Phase']
        # Aggregate NP by phase label over the vertex dimension
        if 'vertex' not in pf.dims or 'Phase' not in eq.data_vars:
            print("DEBUG: 'phase' and 'vertex' missing or no 'Phase' var; eq dims:", pf.dims)
            print("DEBUG: eq coords keys:", list(eq.coords))
            print("DEBUG: phases requested:", phases[:15], ("..." if len(phases) > 15 else ""))
            # Emit empty results
            summary = {
                'tdb_used': str(tdb_path),
                'components': comps,
                'temperature_K': 750.0,
                'pressure_Pa': 101325.0,
                'moles_total': 1.0,
                'composition': {'X(AL)': 0.70, 'X(FE)': 0.30},
                'phases_present': [],
                'claim': 'Two-phase mixture of Al13Fe4 and Al-rich BCC (A2) solid solution',
                'claim_correct': False,
                'notes': {
                    'reason': "No 'phase' or 'vertex' available to resolve phase fractions."
                }
            }
            with open(out_json, 'w') as f:
                json.dump(summary, f, indent=2)
            with open(out_csv, 'w', newline='') as f:
                w = csv.writer(f)
                w.writerow(['phase', 'fraction'])
            print("Database:", tdb_path)
            print("Temperature: 750 K; Pressure: 1 atm; X(AL)=0.70 (X(Fe)=0.30)")
            print("Stable phases (fraction > 1e-3): none")
            print("Claim (Al13Fe4 + BCC_A2 two-phase) correct: False")
            return

        # Sum moles by phase label
        phase_labels = eq['Phase'].values  # dims: N,P,T,X_AL,vertex
        np_vals = pf.values  # same dims
        # Squeeze singleton dims
        phase_labels = np.squeeze(phase_labels)
        np_vals = np.squeeze(np_vals)
        # Ensure 1D over vertex
        # Handle case where there are no vertices (degenerate), then treat as empty
        if np_vals.ndim == 1:
            np_per_vertex = np_vals
            phases_per_vertex = phase_labels
        else:
            # Flatten all but last dim (vertex is last)
            # Our conditions produce shape (vertex,), but keep generality
            flat = np_vals.reshape(-1, np_vals.shape[-1])
            flat_phase = phase_labels.reshape(-1, phase_labels.shape[-1])
            # Take the first row (all rows identical under single state point)
            np_per_vertex = flat[0]
            phases_per_vertex = flat_phase[0]
        phase_fractions = {}
        total_moles = float(np.nansum(np_per_vertex))
        if total_moles <= 0:
            present = []
        else:
            for ph, moles in zip(phases_per_vertex, np_per_vertex):
                if ph is None or ph == '' or not np.isfinite(moles):
                    continue
                phase_fractions[ph] = phase_fractions.get(ph, 0.0) + float(moles)
            present = sorted(phase_fractions.items(), key=lambda x: x[1], reverse=True)

        # Thresholding by fraction of total moles
        thresh = 1e-3
        present = [(p, m/total_moles) for p, m in present if (m/total_moles) > thresh]

        # Summarize claim
        ph_names = {p for p, _ in present}
        has_al13fe4 = any('AL13' in p and 'FE4' in p for p in ph_names) or ('AL13FE4' in ph_names)
        bcc_aliases = {'BCC_A2', 'A2', 'BCC'}
        has_bcc_a2 = any(p in bcc_aliases or 'BCC_A2' in p or (p.endswith('_A2') and 'BCC' in p) for p in ph_names)
        is_two_phase = (len(present) == 2)
        claim_correct = is_two_phase and has_al13fe4 and has_bcc_a2

        summary = {
            'tdb_used': str(tdb_path),
            'components': comps,
            'temperature_K': 750.0,
            'pressure_Pa': 101325.0,
            'moles_total': 1.0,
            'composition': {'X(AL)': 0.70, 'X(FE)': 0.30},
            'phases_present': [{'phase': p, 'fraction': f} for p, f in present],
            'claim': 'Two-phase mixture of Al13Fe4 and Al-rich BCC (A2) solid solution',
            'claim_correct': bool(claim_correct),
            'notes': {
                'aggregation': 'Fractions derived by summing NP over vertices grouped by Phase label.'
            }
        }
        with open(out_json, 'w') as f:
            json.dump(summary, f, indent=2)
        with open(out_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['phase', 'fraction'])
            for p, frac in present:
                w.writerow([p, f"{frac:.6f}"])
        print("Database:", tdb_path)
        print("Temperature: 750 K; Pressure: 1 atm; X(AL)=0.70 (X(Fe)=0.30)")
        if not present:
            print("Stable phases (fraction > 1e-3): none")
        else:
            print("Stable phases (fraction > 1e-3):")
            for p, frac in present:
                print(f"  - {p}: {frac:.6f}")
        print(f"Claim (Al13Fe4 + BCC_A2 two-phase) correct: {claim_correct}")
        return

    tot = pf.sum('phase')
    phase_frac = (pf / tot)

    # Extract fractions for this single state point
    # Set a small threshold for presence
    thresh = 1e-3
    present = []
    for ph in phase_frac['phase'].values.tolist():
        frac = float(phase_frac.sel(phase=ph).values.squeeze())
        if frac > thresh:
            present.append((ph, frac))

    # Sort by fraction descending
    present.sort(key=lambda x: x[1], reverse=True)

    # Summarize claim: two-phase of AL13FE4 and BCC_A2 (Al-rich)
    ph_names = {p for p, _ in present}
    has_al13fe4 = any('AL13' in p and 'FE4' in p for p in ph_names) or ('AL13FE4' in ph_names)
    # Common labels for BCC disordered
    bcc_aliases = {'BCC_A2', 'A2', 'BCC'}
    has_bcc_a2 = any(p in bcc_aliases or 'BCC_A2' in p or (p.endswith('_A2') and 'BCC' in p) for p in ph_names)

    is_two_phase = (len(present) == 2)
    claim_correct = is_two_phase and has_al13fe4 and has_bcc_a2

    summary = {
        'tdb_used': str(tdb_path),
        'components': comps,
        'temperature_K': 750.0,
        'pressure_Pa': 101325.0,
        'moles_total': 1.0,
        'composition': {'X(AL)': 0.70, 'X(FE)': 0.30},
        'phases_present': [{'phase': p, 'fraction': f} for p, f in present],
        'claim': 'Two-phase mixture of Al13Fe4 and Al-rich BCC (A2) solid solution',
        'claim_correct': bool(claim_correct),
        'notes': {
            'has_AL13FE4': bool(has_al13fe4),
            'has_BCC_A2': bool(has_bcc_a2),
            'num_phases_above_threshold': len(present)
        }
    }

    # Write JSON
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)

    # Write CSV of fractions
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['phase', 'fraction'])
        for p, frac in present:
            w.writerow([p, f"{frac:.6f}"])

    # Print concise human-readable summary
    print("Database:", tdb_path)
    print("Temperature: 750 K; Pressure: 1 atm; X(AL)=0.70 (X(Fe)=0.30)")
    if not present:
        print("No phases found above threshold; check database coverage for Fe–Al.")
    else:
        print("Stable phases (fraction > 1e-3):")
        for p, frac in present:
            print(f"  - {p}: {frac:.6f}")
    print(f"Claim (Al13Fe4 + BCC_A2 two-phase) correct: {claim_correct}")


if __name__ == '__main__':
    np.set_printoptions(suppress=True)
    main()
