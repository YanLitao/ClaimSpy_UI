import json
import math
import os
from typing import Dict, List, Tuple

import numpy as np
from pycalphad import Database, equilibrium, variables as v


def normalize(name: str) -> str:
    return ''.join(ch for ch in name.lower() if ch.isalnum())


def find_compound_phase(phases: List[str], target_formula: str, components: Tuple[str, str, str]) -> str:
    """
    Search for a phase name that encodes the target_formula (e.g., 'TI4NB3AL9') allowing
    for permutations of the element order and ignoring separators/underscore/case.
    Returns the matching phase name if found, else ''.
    """
    els = [c.upper() for c in components]
    # Generate all 6 permutations of the element order in the formula string
    perms = [
        target_formula,
        target_formula.replace(els[0], els[1]).replace(els[1], els[0]),
    ]
    # Build all permutations explicitly for safety
    ti, nb, al = els
    formulas = [
        f"{ti}4{nb}3{al}9",
        f"{ti}4{al}9{nb}3",
        f"{nb}3{ti}4{al}9",
        f"{nb}3{al}9{ti}4",
        f"{al}9{ti}4{nb}3",
        f"{al}9{nb}3{ti}4",
    ]
    targets = {normalize(f) for f in formulas}
    for ph in phases:
        nph = normalize(ph)
        if any(t in nph for t in targets):
            return ph
    return ''


def pick_bcc_phase(phases: List[str]) -> str:
    # Prefer canonical BCC_A2 if present
    for ph in phases:
        up = ph.upper()
        if 'BCC' in up and 'A2' in up:
            return ph
    # Fallbacks: names containing 'BETA' or 'TI_BCC'
    for ph in phases:
        up = ph.upper()
        if 'BETA' in up or ('BCC' in up):
            return ph
    return ''


def main():
    # Configuration: local-only TDB directory (from CALPHAD.md guidance)
    tdb_dir = '/Users/delip/play/miniclaimspy/TnETDBDB'
    # Discover available TDBs (local only)
    available_tdbs = [os.path.join(tdb_dir, f) for f in os.listdir(tdb_dir) if f.lower().endswith('.tdb')]
    # Record considered and chosen TDB
    considered = available_tdbs.copy()
    if not available_tdbs:
        raise RuntimeError(f"No TDB files found in {tdb_dir}")

    # For this environment there is typically only COST507R; use it if present
    chosen = ''
    for path in available_tdbs:
        if os.path.basename(path).lower() == 'cost507r.tdb':
            chosen = path
            break
    if not chosen:
        # If not found, just pick the first, but we will validate elements
        chosen = available_tdbs[0]

    # Load database
    dbf = Database(chosen)

    # Elements and phases
    elements = {el.upper() for el in dbf.elements if el != 'VA'}
    phases = sorted(dbf.phases.keys())

    # Required elements
    need = {'TI', 'AL', 'NB'}
    have_all_elements = need.issubset(elements)

    # Phase name discovery
    compound_target = 'TI4NB3AL9'
    compound_phase = find_compound_phase(phases, compound_target, ('TI', 'NB', 'AL'))
    bcc_phase = pick_bcc_phase(phases)

    # Prepare result structure
    result: Dict = {
        'tdb': {
            'chosen': chosen,
            'considered': considered,
        },
        'database_elements': sorted(elements),
        'database_phases_count': len(phases),
        'queries': {
            'compound_name_query': compound_target,
            'compound_phase_found': bool(compound_phase),
            'compound_phase_name': compound_phase,
            'bcc_phase_found': bool(bcc_phase),
            'bcc_phase_name': bcc_phase,
            'elements_ok': have_all_elements,
        },
        'conditions': {
            'T_K': 1273.0,
            'P_Pa': 101325.0,
        },
        'two_phase_field': {
            'exists': False,
            'threshold': 1e-3,
            'composition_range': None,  # to be filled with {'X_TI': [min,max], 'X_AL': [...], 'X_NB': [...]}
            'num_grid_points': 0,
        },
    }

    # Quick fail conditions
    if not have_all_elements:
        print('TDB suitability check: required elements TI/AL/NB are not all present.')
        print(f"Elements in DB: {sorted(elements)}")
        # Write results and exit
        with open('results.json', 'w') as f:
            json.dump(result, f, indent=2)
        # Minimal CSV describing emptiness
        with open('results.csv', 'w') as f:
            f.write('X_AL,X_NB,X_TI,frac_compound,frac_bcc,stable_compound,stable_bcc\n')
        print('Summary: Database unsuitable for Ti–Al–Nb; aborting equilibrium calculation.')
        return

    if not compound_phase:
        print(f"Compound phase '{compound_target}' not found in database phases.")
        with open('results.json', 'w') as f:
            json.dump(result, f, indent=2)
        with open('results.csv', 'w') as f:
            f.write('X_AL,X_NB,X_TI,frac_compound,frac_bcc,stable_compound,stable_bcc\n')
        print('Summary: Target compound absent; cannot form a two-phase field with β-Ti.')
        return

    if not bcc_phase:
        print('No BCC/β-Ti phase identified (e.g., BCC_A2) in this database.')
        with open('results.json', 'w') as f:
            json.dump(result, f, indent=2)
        with open('results.csv', 'w') as f:
            f.write('X_AL,X_NB,X_TI,frac_compound,frac_bcc,stable_compound,stable_bcc\n')
        print('Summary: β-Ti absent; cannot form the queried two-phase field.')
        return

    # Set conditions: 1273 K isothermal section around stoichiometric composition of Ti4Nb3Al9
    x_ti0, x_nb0, x_al0 = 4/16, 3/16, 9/16
    # Explore a modest window around stoichiometry to detect coexistence; bounds clipped to [eps, 1-eps]
    half_width = 0.05
    step = 0.01
    x_al_lo, x_al_hi = max(1e-6, x_al0 - half_width), min(1 - 1e-6, x_al0 + half_width)
    x_nb_lo, x_nb_hi = max(1e-6, x_nb0 - half_width), min(1 - 1e-6, x_nb0 + half_width)

    x_al_grid = np.round(np.arange(x_al_lo, x_al_hi + 1e-12, step), 5)
    x_nb_grid = np.round(np.arange(x_nb_lo, x_nb_hi + 1e-12, step), 5)

    # Broadcasted grid over X(AL) and X(NB); X(TI) implied
    conds = {
        v.T: 1273.0,
        v.P: 101325.0,
        v.N: 1.0,
        v.X('AL'): x_al_grid,
        v.X('NB'): x_nb_grid,
    }

    # Compute equilibrium using all phases (safer to avoid biasing stability)
    comps = ['TI', 'AL', 'NB']
    eq = equilibrium(dbf, comps, phases, conds, verbose=False)

    pf = eq['NP']
    if 'vertex' in pf.dims:
        pf = pf.sum('vertex')
    tot = pf.sum('phase')
    phase_frac = pf / tot

    # Mask for target phases exceeding threshold
    thr = 1e-3
    # Align dims: phase_frac dims typically ('T','P','X_AL','X_NB','phase')
    dims = phase_frac.dims
    # Build boolean arrays for each phase
    phase_names = list(phase_frac['phase'].values)
    try:
        idx_comp = phase_names.index(compound_phase)
        idx_bcc = phase_names.index(bcc_phase)
    except ValueError:
        # In case of name mismatch due to pycalphad phase label normalization
        # try case-insensitive match
        def find_idx(target: str) -> int:
            target_up = target.upper()
            for i, p in enumerate(phase_names):
                if p.upper() == target_up:
                    return i
            raise
        idx_comp = find_idx(compound_phase)
        idx_bcc = find_idx(bcc_phase)

    frac_comp = phase_frac.isel(phase=idx_comp)
    frac_bcc = phase_frac.isel(phase=idx_bcc)

    mask_comp = (frac_comp > thr)
    mask_bcc = (frac_bcc > thr)
    mask_both = mask_comp & mask_bcc

    # Reduce over all non-composition dims except X_AL and X_NB
    reduce_dims = set(mask_both.dims) - {'X_AL', 'X_NB'}
    for d in list(reduce_dims):
        mask_both = mask_both.any(dim=d)
        frac_comp = frac_comp.mean(dim=d)
        frac_bcc = frac_bcc.mean(dim=d)

    # Extract composition points where both phases coexist
    coexist_points: List[Tuple[float, float, float]] = []
    co_fracs: List[Tuple[float, float]] = []
    al_vals = mask_both.coords['X_AL'].values
    nb_vals = mask_both.coords['X_NB'].values

    for i, xa in enumerate(al_vals):
        for j, xn in enumerate(nb_vals):
            if bool(mask_both.values[i, j]):
                xt = 1.0 - float(xa) - float(xn)
                if xt >= 0.0:  # valid composition simplex
                    coexist_points.append((float(xt), float(xa), float(xn)))
                    co_fracs.append((float(frac_comp.values[i, j]), float(frac_bcc.values[i, j])))

    exists = len(coexist_points) > 0
    result['two_phase_field']['exists'] = exists
    result['two_phase_field']['num_grid_points'] = len(coexist_points)

    if exists:
        x_ti_vals = [p[0] for p in coexist_points]
        x_al_vals = [p[1] for p in coexist_points]
        x_nb_vals = [p[2] for p in coexist_points]
        comp_range = {
            'X_TI': [min(x_ti_vals), max(x_ti_vals)],
            'X_AL': [min(x_al_vals), max(x_al_vals)],
            'X_NB': [min(x_nb_vals), max(x_nb_vals)],
        }
        result['two_phase_field']['composition_range'] = comp_range

    # Write CSV of the sampled grid with phase fractions for the two target phases
    with open('results.csv', 'w') as f:
        f.write('X_AL,X_NB,X_TI,frac_compound,frac_bcc,stable_compound,stable_bcc\n')
        k = 0
        for i, xa in enumerate(al_vals):
            for j, xn in enumerate(nb_vals):
                xt = 1.0 - float(xa) - float(xn)
                if xt < 0:
                    continue
                fc = float(frac_comp.values[i, j])
                fb = float(frac_bcc.values[i, j])
                sc = 1 if fc > thr else 0
                sb = 1 if fb > thr else 0
                f.write(f"{xa:.5f},{xn:.5f},{xt:.5f},{fc:.6g},{fb:.6g},{sc},{sb}\n")
                k += 1

    # Write JSON summary
    with open('results.json', 'w') as f:
        json.dump(result, f, indent=2)

    # Print concise human summary
    print('--- CALPHAD Summary (Ti–Al–Nb at 1273 K) ---')
    print(f"TDB chosen: {chosen}")
    print(f"Elements in DB: {sorted(elements)}")
    print(f"Compound phase present: {bool(compound_phase)} ({compound_phase})")
    print(f"β-Ti (BCC) phase present: {bool(bcc_phase)} ({bcc_phase})")
    if exists:
        cr = result['two_phase_field']['composition_range']
        print('Two-phase field (compound + β-Ti): PRESENT within sampled window.')
        print(f"Approx. composition bounds (mole fraction): ")
        print(f"  X(TI): {cr['X_TI'][0]:.4f} – {cr['X_TI'][1]:.4f}")
        print(f"  X(AL): {cr['X_AL'][0]:.4f} – {cr['X_AL'][1]:.4f}")
        print(f"  X(NB): {cr['X_NB'][0]:.4f} – {cr['X_NB'][1]:.4f}")
        print(f"Grid points with coexistence: {result['two_phase_field']['num_grid_points']}")
    else:
        print('Two-phase field (compound + β-Ti): NOT FOUND within the sampled window around stoichiometry.')


if __name__ == '__main__':
    main()

