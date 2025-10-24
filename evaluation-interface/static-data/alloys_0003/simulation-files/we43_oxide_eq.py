import json, csv, sys, os
from pathlib import Path
import numpy as np
from pycalphad import Database, equilibrium, variables as v

# ---- Configuration (single-run script per CALPHAD.md guidance) ----
# Paths
CWD = Path.cwd()
TDB_PATH = Path('/Users/delip/play/miniclaimspy/TnETDBDB/cost507R.TDB')
assert TDB_PATH.exists(), f"TDB not found: {TDB_PATH}"

# Alloy: WE43 with 1.0 wt% O; balance Mg
# Composition in weight percent (normalized to 100 wt%)
wt_comp = {
    'MG': 91.5,
    'Y': 4.0,
    'ND': 3.0,
    'ZR': 0.5,
    'O': 1.0,
}

# Atomic weights (g/mol)
MASS = {
    'MG': 24.305,
    'Y': 88.90584,
    'ND': 144.242,
    'ZR': 91.224,
    'O': 15.999,
}

# Convert weight% to mole fractions (exclude VA)
# Normalize to 1.0 total
moles = {el: wt_comp[el] / MASS[el] for el in wt_comp}
tot_moles = sum(moles.values())
x = {el: moles[el] / tot_moles for el in moles}
# Enforce small numerical tolerance and renormalize
neg = [k for k, v_ in x.items() if v_ < 0]
assert not neg, f"Negative mole fractions computed for {neg}"
sumx = sum(x.values())
for k in x:
    x[k] /= sumx

# Components (include VA for solution models)
components = ['MG', 'Y', 'ND', 'ZR', 'O', 'VA']

# Conditions: N−1 mole fraction constraints for N substitutional components (exclude VA)
conds = {
    v.T: np.array([298.0, 900.0]),
    v.P: 101325.0,
    v.N: 1.0,
    v.X('MG'): float(x['MG']),
    v.X('Y'): float(x['Y']),
    v.X('ND'): float(x['ND']),
    v.X('ZR'): float(x['ZR']),
    # X('O') is implied to satisfy sum=1 among substitutional elements
}

# Load database and select phases containing at least one of our components
print(f"Loading database: {TDB_PATH}")
dbf = Database(str(TDB_PATH))

sel_comps = set(['MG', 'Y', 'ND', 'ZR', 'O'])
phases = []
for ph in sorted(dbf.phases.keys()):
    const = dbf.phases[ph].constituents
    # Flatten constituents to element symbols
    flat = set()
    for subl in const:
        for c in subl:
            flat.add(str(c))
    if flat & sel_comps:
        phases.append(ph)

assert phases, "No phases selected; check database contents."
print(f"Selected {len(phases)} phases relevant to components {sorted(sel_comps)}")

# Run equilibrium
print("Computing equilibrium at T=298 K and 900 K ...")
eq = equilibrium(dbf, components, phases, conds, verbose=False)

"""
In pycalphad, when no explicit 'phase' dimension exists, the equilibrium is
represented as a convex combination over 'vertex' entries. Each vertex has a
phase label in eq['Phase'] and a corresponding amount in eq['NP'].
Aggregate NP by Phase across vertices to get per-phase moles.
"""

# Identify oxide phases by presence of 'O' in their constituents
oxide_phases = set()
for ph in phases:
    const = dbf.phases[ph].constituents
    if any('O' in {str(c) for c in subl} for subl in const):
        oxide_phases.add(ph)

# Prepare results for each temperature
temps = list(eq.T.values)
res_rows = []
summary = {
    'composition_wt%': wt_comp,
    'composition_mole_fraction': {k: float(v) for k, v in x.items()},
    'temperatures_K': temps,
    'pressure_Pa': 101325.0,
    'tdb': str(TDB_PATH),
    'notes': 'Phase fractions are mole-based. Oxides identified via constituent O in phase model. Aggregation uses NP-by-vertex grouped by Phase.'
}

stable_by_T = {}
for iT, T in enumerate(temps):
    # Slice NP and Phase at this T and current composition
    np_sel = eq['NP'].isel(T=iT, P=0, N=0, X_MG=0, X_Y=0, X_ND=0, X_ZR=0)
    ph_sel = eq['Phase'].isel(T=iT, P=0, N=0, X_MG=0, X_Y=0, X_ND=0, X_ZR=0)

    # Aggregate NP by Phase label
    per_phase = {}
    total_moles = float(np_sel.sum('vertex').values)
    # Guard against zero total (should not happen)
    if not np.isfinite(total_moles) or total_moles <= 0:
        total_moles = 1.0

    for iv in range(np_sel.sizes['vertex']):
        fmol = float(np_sel.isel(vertex=iv).values)
        if fmol <= 0:
            continue
        ph_name = str(ph_sel.isel(vertex=iv).values)
        # Some vertices may be empty (phase name 'None'); skip
        if not ph_name or ph_name.lower() == 'none':
            continue
        per_phase[ph_name] = per_phase.get(ph_name, 0.0) + fmol

    # Convert to fractions
    stable = []
    stable_oxides = []
    for ph_name, nmol in per_phase.items():
        frac = nmol / total_moles
        if frac > 1e-8:
            is_ox = ph_name in oxide_phases
            res_rows.append({
                'T_K': float(T),
                'phase': ph_name,
                'mole_fraction': float(frac),
                'is_oxide': bool(is_ox)
            })
            stable.append((ph_name, float(frac)))
            if is_ox:
                stable_oxides.append((ph_name, float(frac)))

    stable.sort(key=lambda t: t[1], reverse=True)
    stable_oxides.sort(key=lambda t: t[1], reverse=True)
    stable_by_T[str(float(T))] = {
        'stable_phases': stable,
        'oxide_phases': stable_oxides,
        'total_oxide_fraction': float(sum(f for _, f in stable_oxides))
    }

summary['stable_by_T'] = stable_by_T

# Write outputs
json_path = CWD / 'we43_oxide_results.json'
csv_path = CWD / 'we43_oxide_results.csv'
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

with open(csv_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['T_K', 'phase', 'mole_fraction', 'is_oxide'])
    w.writeheader()
    for row in res_rows:
        w.writerow(row)

# Print concise human summary
print("=== CALPHAD Summary ===")
print(f"Alloy (wt%): {wt_comp}")
print("Mole fractions:")
for k in ['MG','Y','ND','ZR','O']:
    print(f"  X({k}) = {x[k]:.6f}")
for T in temps:
    info = stable_by_T[str(float(T))]
    top = ", ".join([f"{p}:{f:.3f}" for p, f in info['stable_phases'][:6]])
    ox = ", ".join([f"{p}:{f:.3f}" for p, f in info['oxide_phases']])
    print(f"T = {T:.1f} K")
    print(f"  Top phases: {top if top else 'None'}")
    print(f"  Oxide phases: {ox if ox else 'None'}")
    print(f"  Total oxide fraction (mole): {info['total_oxide_fraction']:.4f}")

# ---- Metal-only equilibrium (WE43 baseline without O) ----
# Rationale: COST507R lacks explicit oxide phases; oxygen cannot be accommodated in metallic phases
# so we also report the metal-only equilibrium for base WE43 (Mg-4Y-3Nd-0.5Zr wt%).

wt_comp_metal = {
    'MG': 92.5,
    'Y': 4.0,
    'ND': 3.0,
    'ZR': 0.5,
}
mm_metal = {el: wt_comp_metal[el] / MASS[el] for el in wt_comp_metal}
tot_m_metal = sum(mm_metal.values())
x_metal = {el: mm_metal[el] / tot_m_metal for el in mm_metal}

components_metal = ['MG', 'Y', 'ND', 'ZR', 'VA']
conds_metal = {
    v.T: np.array([298.0, 900.0]),
    v.P: 101325.0,
    v.N: 1.0,
    # For 4 substitutional components, set 3 constraints; X(ZR) implied
    v.X('MG'): float(x_metal['MG']),
    v.X('Y'): float(x_metal['Y']),
    v.X('ND'): float(x_metal['ND']),
}

# Select phases that contain any of these metals
phases_metal = []
for ph in sorted(dbf.phases.keys()):
    const = dbf.phases[ph].constituents
    flat = set()
    for subl in const:
        for c in subl:
            flat.add(str(c))
    if flat & set(['MG','Y','ND','ZR']):
        phases_metal.append(ph)

eq_metal = equilibrium(dbf, components_metal, phases_metal, conds_metal, verbose=False)

def per_phase_fractions_from_eq(eqds, comp_keys):
    temps_local = list(eqds.T.values)
    out = {}
    for iT, T in enumerate(temps_local):
        np_sel = eqds['NP'].isel(T=iT, P=0, N=0, **{f'X_{k}':0 for k in comp_keys})
        ph_sel = eqds['Phase'].isel(T=iT, P=0, N=0, **{f'X_{k}':0 for k in comp_keys})
        total = float(np_sel.sum('vertex').values)
        if not np.isfinite(total) or total <= 0:
            phases_ = []
        else:
            accum = {}
            for iv in range(np_sel.sizes['vertex']):
                nmol = float(np_sel.isel(vertex=iv).values)
                if nmol <= 0: continue
                pname = str(ph_sel.isel(vertex=iv).values)
                if not pname or pname.lower()=='none':
                    continue
                accum[pname] = accum.get(pname, 0.0) + nmol
            phases_ = sorted([(k, v/total) for k, v in accum.items()], key=lambda t: t[1], reverse=True)
        out[str(float(T))] = phases_
    return temps_local, out

# Build comp_keys present in eq_metal dims (exclude 'ZR' if not explicit)
dim_keys = [d for d in eq_metal.dims if d.startswith('X_')]
comp_keys = [d.split('_',1)[1] for d in dim_keys]
temps_metal, stable_metal_by_T = per_phase_fractions_from_eq(eq_metal, comp_keys)

# Persist metal-only results
json_path_metal = CWD / 'we43_metal_results.json'
csv_path_metal = CWD / 'we43_metal_results.csv'
metal_summary = {
    'composition_wt%': wt_comp_metal,
    'composition_mole_fraction': {k: float(v) for k, v in x_metal.items()},
    'temperatures_K': temps_metal,
    'pressure_Pa': 101325.0,
    'tdb': str(TDB_PATH),
    'stable_by_T': stable_metal_by_T,
}
with open(json_path_metal, 'w') as f:
    json.dump(metal_summary, f, indent=2)

with open(csv_path_metal, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['T_K','phase','mole_fraction'])
    for Tstr, plist in stable_metal_by_T.items():
        for ph, frac in plist:
            w.writerow([float(Tstr), ph, float(frac)])

# ---- Stoichiometric oxide estimate (upper-bound) ----
# Allocate 1.0 g O to oxides in order of expected oxygen affinity: Y -> ND -> ZR -> MG
O_avail = 1.0  # g per 100 g alloy
oxide_seq = [
    ('Y', 'Y2O3', 2*MASS['Y'] + 3*MASS['O'], 3*MASS['O']/(2*MASS['Y'])),  # (element, oxide_name, M_oxide, mO_per_g_element)
    ('ND','Nd2O3', 2*MASS['ND'] + 3*MASS['O'], 3*MASS['O']/(2*MASS['ND'])),
    ('ZR','ZrO2', MASS['ZR'] + 2*MASS['O'], 2*MASS['O']/MASS['ZR']),
    ('MG','MgO', MASS['MG'] + MASS['O'], MASS['O']/MASS['MG']),
]
available = {'Y': 4.0, 'ND': 3.0, 'ZR': 0.5, 'MG': 91.5}  # g available in alloy when O=1.0 wt%
formed = []  # list of (oxide, mass_g, element_consumed_g, O_used_g)
for el, ox, Mox, o_per_el in oxide_seq:
    if O_avail <= 0: break
    # Max element that can be oxidized by remaining oxygen
    el_needed_per_O = 1.0 / o_per_el  # g el per g O
    el_needed = min(available[el], O_avail * el_needed_per_O)
    if el_needed <= 0: continue
    O_used = el_needed * o_per_el
    mass_oxide = el_needed + O_used
    available[el] -= el_needed
    O_avail -= O_used
    formed.append((ox, mass_oxide, el_needed, O_used))

stoich = {
    'oxides': [{'oxide': ox, 'mass_g_per_100g': mass, 'element_consumed_g': elc, 'O_used_g': Ou}
               for (ox, mass, elc, Ou) in formed],
    'total_oxide_mass_g_per_100g': float(sum(m for (_, m, _, _) in formed)),
    'O_unreacted_g_per_100g': float(O_avail),
}

# Merge stoichiometric estimate into the oxygen-run JSON
summary['stoichiometric_oxide_estimate'] = stoich
with open(json_path, 'w') as f:
    json.dump(summary, f, indent=2)

print("\n=== Metal-only Equilibrium (no O) ===")
print(f"Alloy (wt%): {wt_comp_metal}")
for k in ['MG','Y','ND','ZR']:
    print(f"  X({k}) = {x_metal[k]:.6f}")
for T in temps_metal:
    plist = stable_metal_by_T[str(float(T))]
    top = ", ".join([f"{p}:{f:.3f}" for p, f in plist[:8]])
    print(f"T = {T:.1f} K")
    print(f"  Top phases: {top if top else 'None'}")

print("\n=== Stoichiometric Oxide Estimate (with 1 wt% O) ===")
for ox in stoich['oxides']:
    print(f"  {ox['oxide']}: {ox['mass_g_per_100g']:.3f} g per 100 g alloy")
print(f"  Total oxides: {stoich['total_oxide_mass_g_per_100g']:.3f} g per 100 g alloy")
if stoich['O_unreacted_g_per_100g'] > 1e-6:
    print(f"  Unreacted O: {stoich['O_unreacted_g_per_100g']:.3f} g per 100 g")
