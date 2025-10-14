#!/usr/bin/env python3
"""
Compute equilibrium phases and fractions for Al–Mg–Zn at 298 K
using the local COST507 (cost507R.TDB) database in TnETDBDB.

Outputs:
- results.json: structured results with phases and mole fractions for each composition
- results.csv: flat table summary (composition, phase, mole_fraction)
- summary.txt: concise human-readable summary

Notes:
- Uses mole-fraction phase fractions based on NP (moles of phase). Do not use 'W' for constraints.
- Detects presence of tau-like phases by name (e.g., contains 'TAU') and reports.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from pycalphad import Database, equilibrium, variables as v


def phase_has_any(dbf: Database, phase: str, comps: List[str]) -> bool:
    pobj = dbf.phases[phase]
    const = pobj.constituents
    species = {str(sp) for subl in const for sp in subl}
    return bool(set(comps) & species)


def get_relevant_phases(dbf: Database, comps: List[str]) -> List[str]:
    phases = []
    for p in dbf.phases:
        if phase_has_any(dbf, p, comps):
            phases.append(p)
    phases = sorted(set(phases))
    return phases


def compute_phase_fractions(eq) -> Tuple[pd.Series, pd.Series]:
    """Return DataFrame of phase mole fractions and Series of totals.

    The eq Dataset may include a 'vertex' dimension; sum it first.
    """
    pf = eq['NP']
    dims = list(pf.dims)
    # We expect NP to have a 'vertex' dimension with per-vertex moles; sum over all other dims
    if 'vertex' not in dims:
        # No vertex dimension; treat as single-point, single-vertex
        total_moles = float(np.atleast_1d(pf.values).sum())
        s = pd.Series({'UNKNOWN': 1.0 if total_moles > 0 else 0.0})
        return s, pd.Series({'total_moles': total_moles})

    # Sum over all dims except 'vertex' to get per-vertex moles
    sum_dims = [d for d in dims if d != 'vertex']
    pf_vertex = pf.sum(sum_dims)

    # Extract per-vertex phase names from eq['Phase'] by selecting the single condition indices
    ph_da = eq['Phase']
    # Reduce all non-vertex dims to the first (and only) index
    for d in ph_da.dims:
        if d != 'vertex':
            ph_da = ph_da.isel({d: 0})
    # Now ph_da should be 1D over 'vertex'
    phase_names = []
    for val in np.atleast_1d(ph_da.values):
        if isinstance(val, (bytes, bytearray)):
            phase_names.append(val.decode())
        else:
            phase_names.append(str(val))

    # Aggregate moles by phase name across vertices
    mol_by_phase: Dict[str, float] = {}
    for i, pname in enumerate(phase_names):
        moles_i = float(np.atleast_1d(pf_vertex.sel(vertex=i).values).squeeze())
        mol_by_phase[pname] = mol_by_phase.get(pname, 0.0) + moles_i

    total_moles = sum(mol_by_phase.values())
    # Convert to fractions; avoid divide-by-zero
    if total_moles > 0:
        frac_by_phase = {p: (m / total_moles) for p, m in mol_by_phase.items()}
    else:
        frac_by_phase = {p: 0.0 for p in mol_by_phase}

    s = pd.Series(frac_by_phase)
    s = s.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return s.sort_values(ascending=False), pd.Series({'total_moles': total_moles})


def main():
    outdir = Path('.')
    tdb_path = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
    assert tdb_path.exists(), f"TDB not found: {tdb_path}"

    # System definition
    comps_main = ['AL', 'MG', 'ZN']
    comps = comps_main + ['VA']  # include vacancies for substitutional models
    T = 298.0  # K
    P = 101325  # Pa

    # Load database exactly once; avoid repeated parsing
    dbf = Database(str(tdb_path))

    # Select relevant phases programmatically by constituents
    phases = get_relevant_phases(dbf, comps_main)
    assert phases, "No phases selected; did you build the phases list?"

    # Detect if a tau-like phase exists in the database nomenclature
    tau_like_phases = sorted([p for p in phases if 'TAU' in p.upper()])

    # Define target compositions (mole fractions)
    # Specify N−1 mole fraction constraints for N components
    compositions = [
        {"label": "Al-4Mg-2Zn at.%", "X_AL": 0.94, "X_MG": 0.04},
        {"label": "Al-7Mg-3Zn at.%", "X_AL": 0.90, "X_MG": 0.07},
        {"label": "Al-2Mg-1Zn at.%", "X_AL": 0.97, "X_MG": 0.02},
    ]

    results: Dict[str, Dict] = {
        'tdb': str(tdb_path),
        'temperature_K': T,
        'pressure_Pa': P,
        'tau_phase_names_in_db': tau_like_phases,
        'entries': [],
    }

    rows: List[Dict] = []

    for comp in compositions:
        label = comp['label']
        x_al = comp['X_AL']
        x_mg = comp['X_MG']
        # Validate N−1 constraints
        if any(v < 0 or v > 1 for v in (x_al, x_mg)) or (x_al + x_mg > 1.0 + 1e-12):
            raise ValueError(f"Invalid composition for {label}: X_AL={x_al}, X_MG={x_mg}")

        conds = {v.T: T, v.P: P, v.N: 1.0, v.X('AL'): x_al, v.X('MG'): x_mg}

        # Compute equilibrium; use higher pdens for robustness on complex phase sets
        try:
            eq = equilibrium(dbf, comps, phases, conds, calc_opts={'pdens': 1500})
        except Exception:
            # fallback to a simpler phase set; expand if needed
            simple = [p for p in phases if p in {'LIQUID', 'FCC_A1', 'HCP_A3', 'HCP_ZN', 'BCC_A2'}]
            if not simple:
                simple = None
            eq = equilibrium(dbf, comps, simple, conds)

        # Minimal debug info about the equilibrium dataset structure (dims and coords)
        dbg = {
            'label': label,
            'NP_dims': tuple(eq['NP'].dims),
            'coords': {k: tuple(v.dims) for k, v in eq.coords.items() if hasattr(v, 'dims')},
            'data_vars': list(eq.data_vars.keys()),
        }
        # Append to a debug file (overwrite on first entry)
        debug_path = outdir / 'debug_eq_structure.jsonl'
        with debug_path.open('a') as f:
            f.write(json.dumps(dbg) + "\n")

        phase_frac, totals = compute_phase_fractions(eq)

        # Threshold for reporting presence
        PRES_TOL = 1e-3
        present = phase_frac[phase_frac > PRES_TOL].sort_values(ascending=False)

        # Primary phase = highest fraction
        primary_phase = present.idxmax() if not present.empty else None
        primary_frac = float(present.iloc[0]) if not present.empty else 0.0

        # Tau detection among present phases
        tau_present = any('TAU' in p.upper() for p in present.index)

        # FCC confirmation
        fcc_present = any(p.upper() == 'FCC_A1' for p in present.index)
        fcc_is_primary = (primary_phase is not None and primary_phase.upper() == 'FCC_A1')

        entry = {
            'label': label,
            'X_AL': x_al,
            'X_MG': x_mg,
            'X_ZN': 1.0 - x_al - x_mg,
            'phases': [
                {'name': p, 'mole_fraction': float(a)} for p, a in present.items()
            ],
            'primary_phase': primary_phase,
            'primary_mole_fraction': primary_frac,
            'tau_present': tau_present,
            'fcc_present': fcc_present,
            'fcc_is_primary': fcc_is_primary,
        }
        results['entries'].append(entry)

        for p, a in present.items():
            rows.append({
                'label': label,
                'phase': p,
                'mole_fraction': float(a),
            })

    # Write JSON and CSV
    (outdir / 'results.json').write_text(json.dumps(results, indent=2))
    pd.DataFrame(rows).to_csv(outdir / 'results.csv', index=False)

    # Concise human summary
    lines: List[str] = []
    lines.append(f"Using TDB: {tdb_path}")
    if tau_like_phases:
        lines.append(f"Tau-like phase names present in DB: {', '.join(tau_like_phases)}")
    else:
        lines.append("No TAU*-named phase is defined in the selected database.")
    lines.append(f"Conditions: T = {T:.1f} K; P = {P} Pa")
    lines.append("")
    for e in results['entries']:
        lines.append(f"Composition {e['label']} (X_AL={e['X_AL']:.3f}, X_MG={e['X_MG']:.3f}, X_ZN={e['X_ZN']:.3f})")
        if e['phases']:
            for ph in e['phases']:
                lines.append(f"  - {ph['name']}: mole fraction {ph['mole_fraction']:.5f}")
        else:
            lines.append("  - No phases above 1e-3 fraction")
        lines.append(f"  Primary phase: {e['primary_phase']} (fraction {e['primary_mole_fraction']:.5f})")
        lines.append(f"  FCC_A1 present: {e['fcc_present']}; FCC_A1 is primary: {e['fcc_is_primary']}")
        lines.append(f"  Tau present among stable phases: {e['tau_present']}")
        lines.append("")
    (outdir / 'summary.txt').write_text("\n".join(lines))

    # Also print succinct summary to stdout
    print("\n".join(lines))


if __name__ == '__main__':
    main()
