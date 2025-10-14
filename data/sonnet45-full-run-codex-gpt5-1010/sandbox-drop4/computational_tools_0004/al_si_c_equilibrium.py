import json
import sys
import warnings
from pathlib import Path

# Keep output concise
warnings.filterwarnings('ignore')

# Load .env (for MP_API_KEY if needed by future workflows)
try:
    from dotenv import load_dotenv
    load_dotenv(Path('.env'))
except Exception:
    pass

from pycalphad import Database, variables as v, equilibrium
import numpy as np


def species_names_in_phase(dbf: Database, phase: str):
    names = set()
    for subl in dbf.phases[phase].constituents:
        for sp in subl:
            try:
                names.add(sp.name)
            except Exception:
                names.add(str(sp))
    return names


def union_species_names(dbf: Database):
    u = set()
    for ph in dbf.phases:
        u |= species_names_in_phase(dbf, ph)
    return u


def main():
    CWD = Path.cwd()
    TDB_DIR = Path('/Users/delip/play/miniclaimspy/TnETDBDB').resolve()
    if not TDB_DIR.exists():
        print(f"ERROR: TDB directory not found: {TDB_DIR}")
        sys.exit(2)

    # Target system: Al-Si-C
    components = ['AL', 'SI', 'C', 'VA']
    base_comps = ['AL', 'SI', 'C']
    conds = {
        v.T: np.array([300.0, 400.0, 500.0]),
        v.P: 101325.0,
        v.N: 1.0,
        v.X('AL'): 0.30,
        v.X('SI'): 0.55,
    }

    # Find a suitable TDB (prefer COST507R)
    all_tdbs = sorted([p for p in TDB_DIR.glob('**/*.TDB')] + [p for p in TDB_DIR.glob('**/*.tdb')])
    if not all_tdbs:
        print(f"ERROR: No TDB files found under {TDB_DIR}")
        sys.exit(2)

    preferred_names = {'cost507R.TDB', 'COST507R.TDB', 'COST507R.tdb', 'cost507r.tdb'}
    preferred = [p for p in all_tdbs if p.name in preferred_names]
    order = preferred + [p for p in all_tdbs if p not in preferred]

    chosen_db = None
    chosen_tdb_path = None
    for tdb_path in order:
        try:
            db = Database(str(tdb_path))
            names = union_species_names(db)
            if set(base_comps).issubset(names):
                chosen_db = db
                chosen_tdb_path = tdb_path
                break
        except Exception:
            continue

    if chosen_db is None:
        print("ERROR: No suitable TDB database found that contains AL, SI, and C.")
        print("Searched in:", TDB_DIR)
        sys.exit(2)

    # Select relevant phases: those that include at least one of AL, SI, C
    def phase_has_any(dbf: Database, phase: str, comps):
        return bool(species_names_in_phase(dbf, phase) & set(comps))

    phases = sorted([p for p in chosen_db.phases if phase_has_any(chosen_db, p, base_comps)])
    assert phases, "No phases selected; did you build the phases list?"

    # Compute equilibrium
    try:
        eq = equilibrium(chosen_db, components, phases, conds, output='GM')
    except Exception as e:
        print("ERROR: Equilibrium calculation failed:", e)
        sys.exit(2)

    # Phase fractions
    pf = eq['NP']  # dims: ('N','P','T','X_AL','X_SI','vertex')
    phase_labels = eq['Phase']  # same dims as pf; string labels per vertex

    Ts = [float(t) for t in eq['T'].values]
    stable_threshold = 1e-6

    # Carbide identification by phase name (robust for common carbides)
    def is_carbide_name(name: str) -> bool:
        n = name.upper()
        if n in {'SIC', 'AL4C3'}:
            return True
        # Common ternaries like AL8SIC7, AL4SI2C5, AL2SIC
        return ('C' in n and ('AL' in n or 'SI' in n))

    carbide_flags = {ph: is_carbide_name(ph) for ph in phases}

    results = {
        'tdb_path': str(chosen_tdb_path),
        'components': components,
        'phases_considered': phases,
        'temperatures_K': Ts,
        'pressure_Pa': float(conds[v.P]),
        'composition': {
            'X_AL': float(conds[v.X('AL')]),
            'X_SI': float(conds[v.X('SI')]),
            'X_C': 1.0 - float(conds[v.X('AL')]) - float(conds[v.X('SI')])
        },
        'phase_fractions': {},
    }

    # Build fractions by grouping NP across vertices with the same Phase label
    for iT, Tval in enumerate(Ts):
        # slice at this temperature
        pf_T = pf.isel(T=iT).squeeze()
        ph_T = phase_labels.isel(T=iT).squeeze()
        # Flatten over remaining dims except 'vertex'
        # Shapes expected: (vertex,)
        np_vals = np.atleast_1d(pf_T.values).reshape(-1)
        labels = np.array([str(x) for x in np.atleast_1d(ph_T.values).reshape(-1)])
        # Aggregate by label
        by_phase = {}
        for val, lbl in zip(np_vals, labels):
            if not lbl or lbl.strip() == '':
                continue
            by_phase[lbl] = by_phase.get(lbl, 0.0) + float(val)
        # Normalize to total
        total = sum(by_phase.values())
        fracs = {ph: (x / total if total > 0 else 0.0) for ph, x in by_phase.items() if (x / total if total > 0 else 0.0) > stable_threshold}
        fracs_sorted = dict(sorted(fracs.items(), key=lambda kv: -kv[1]))
        results['phase_fractions'][str(int(Tval))] = fracs_sorted

    # Emit files
    json_path = CWD / 'al_si_c_results.json'
    csv_path = CWD / 'al_si_c_results.csv'
    summary_path = CWD / 'al_si_c_summary.txt'

    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)

    with open(csv_path, 'w') as f:
        f.write('T_K,phase,mole_fraction,is_carbide\n')
        for Tkey in results['phase_fractions']:
            for ph, x in results['phase_fractions'][Tkey].items():
                f.write(f"{Tkey},{ph},{x:.6f},{'yes' if carbide_flags.get(ph, False) else 'no'}\n")

    # Concise printed summary
    lines = []
    lines.append(f"Using TDB: {chosen_tdb_path}")
    lines.append(f"Components: {components}; Phases considered: {len(phases)}")
    lines.append(f"Composition (at.%): Al 30, Si 55, C 15; P = 101325 Pa")

    all_present = set()
    for Tkey, fr in results['phase_fractions'].items():
        all_present |= set(fr.keys())
    carbides_present = [ph for ph in sorted(all_present) if carbide_flags.get(ph, False)]
    lines.append("Carbide phases present at some T: " + (", ".join(carbides_present) if carbides_present else "None"))

    for Tkey in sorted(results['phase_fractions'].keys(), key=lambda s: int(s)):
        lines.append(f"\nT = {Tkey} K: stable phases (mole fraction)")
        for ph, x in results['phase_fractions'][Tkey].items():
            tag = ' [carbide]' if carbide_flags.get(ph, False) else ''
            lines.append(f"  - {ph}: {x:.6f}{tag}")

    # Dominant carbide info
    carbide_avg = {}
    for ph in carbides_present:
        tot = 0.0
        for Tkey in results['phase_fractions']:
            tot += results['phase_fractions'][Tkey].get(ph, 0.0)
        carbide_avg[ph] = tot / max(1, len(results['phase_fractions']))
    if carbide_avg:
        dom = max(carbide_avg.items(), key=lambda kv: kv[1])[0]
        first_T = None
        for Tkey in sorted(results['phase_fractions'].keys(), key=lambda s: int(s)):
            if results['phase_fractions'][Tkey].get(dom, 0.0) > stable_threshold:
                first_T = int(Tkey)
                break
        lines.append(f"\nDominant carbide: {dom}; first appears by {first_T} K (threshold {stable_threshold}).")

    summary = "\n".join(lines)
    print(summary)
    with open(summary_path, 'w') as f:
        f.write(summary + "\n")


if __name__ == '__main__':
    main()
