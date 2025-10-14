import json
import os
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

from pycalphad import Database, equilibrium, variables as v

# Fixed settings per instructions
T = 973.0  # K
P = 101325  # Pa
X_AL_list = [0.50, 0.55, 0.60, 0.65, 0.70]
NP_TOL = 1e-6

# Choose TDB from the prescribed directory
TDB_DIR = "/Users/delip/play/miniclaimspy/TnETDBDB"
TDB_PATH = os.path.join(TDB_DIR, "cost507R.TDB")
TDB_README = os.path.join(TDB_DIR, "cost507R_README.md")

if not os.path.isfile(TDB_PATH):
    print(f"ERROR: Required TDB not found: {TDB_PATH}")
    sys.exit(2)

# Load database
try:
    dbf = Database(TDB_PATH)
except Exception as e:
    print(f"ERROR: Failed to load TDB {TDB_PATH}: {e}")
    sys.exit(2)

# Components (always include VA for substitutional models)
components = ['AL', 'FE', 'VA']

# Build a relevant phase list based on constituents containing AL or FE
relevant_phases = []
for pname, pobj in dbf.phases.items():
    species = {str(sp) for subl in pobj.constituents for sp in subl}
    if ('AL' in species) or ('FE' in species):
        relevant_phases.append(pname)

# Simple sanity check
if not relevant_phases:
    print("ERROR: No relevant phases found containing AL or FE in the database.")
    sys.exit(2)

# Conditions: specify N-1 mole fraction constraints for substitutional components
conditions = {v.T: T, v.P: P, v.N: 1.0, v.X('AL'): X_AL_list}

# Try broad phase set; if that fails, fallback to a simpler plausible set
fallback_candidates = [p for p in ['LIQUID','FCC_A1','BCC_A2','B2_BCC','B2','B32'] if p in dbf.phases]

calc_opts = {'pdens': 1500}

try:
    eq = equilibrium(dbf, components, relevant_phases, conditions, calc_opts=calc_opts)
except Exception as e:
    if fallback_candidates:
        try:
            eq = equilibrium(dbf, components, fallback_candidates, conditions)
            relevant_phases = fallback_candidates[:]  # record that we used fallback
        except Exception as e2:
            print("ERROR: Equilibrium failed for both full and fallback phase sets.")
            print(f"Primary error: {e}")
            print(f"Fallback error: {e2}")
            sys.exit(2)
    else:
        print("ERROR: Equilibrium failed for full phase set, and no fallback phases available.")
        print(f"Primary error: {e}")
        sys.exit(2)

# Phase fractions per composition
def _decode(name):
    return name.decode() if isinstance(name, (bytes, bytearray)) else str(name)

# Build results per composition using a robust grouping by Phase
results = []
rows = []

def detect_al2fe_like(phname: str) -> bool:
    n = phname.upper()
    targets = [
        'AL2FE','FEAL2',  # FeAl2
        'AL5FE2','FE2AL5',
        'AL13FE4','FE4AL13',
        'AL3FE','FEAL3',
        'ALFE'  # FeAl (B2)
    ]
    return any(t == n or t in n for t in targets)

al2fe_like_present_by_x = {}

# Extract data using groupby by Phase then normalizing per composition
for i, x_al in enumerate(X_AL_list):
    # Select this composition (and collapse singleton dims)
    sel = eq.sel(X_AL=x_al).squeeze()
    # Group NP by Phase and sum over any residual vertices
    grouped = sel['NP'].groupby(sel['Phase']).sum()
    if 'vertex' in grouped.dims:
        grouped = grouped.sum('vertex')
    # Also sum over any leftover dims except the group key 'Phase'
    for d in list(grouped.dims):
        if d not in ('Phase',):
            grouped = grouped.sum(d)
    # Construct fractions
    phases_here = [ _decode(p) for p in list(grouped['Phase'].values) ]
    amounts = grouped.values.astype(float)
    total = float(np.sum(amounts))
    present = []
    al2fe_like_here = False
    if total > 0:
        for pname, amt in zip(phases_here, amounts):
            frac = amt/total
            if frac > 1e-3:
                present.append({"phase": pname, "mole_fraction": float(frac)})
                if detect_al2fe_like(pname):
                    al2fe_like_here = True
    present_sorted = sorted(present, key=lambda d: d['mole_fraction'], reverse=True)
    results.append({
        "X_AL": float(x_al),
        "T_K": float(T),
        "P_Pa": float(P),
        "phases": present_sorted
    })
    al2fe_like_present_by_x[float(x_al)] = al2fe_like_here
    for entry in present_sorted:
        rows.append(OrderedDict(
            X_AL=float(x_al), T_K=float(T), P_Pa=float(P),
            phase=entry['phase'], mole_fraction=entry['mole_fraction']
        ))

# Write outputs in the current working directory
out_json = 'fe_al_eq_973K_results.json'
out_csv = 'fe_al_eq_973K_summary.csv'
provenance = {
    "tdb_path": TDB_PATH,
    "tdb_readme": TDB_README,
    "phase_count_used": len(relevant_phases),
    "phases_used": relevant_phases,
}
with open(out_json, 'w') as f:
    json.dump({
        "system": {"elements": components},
        "conditions": {"T_K": T, "P_Pa": P, "X_AL": X_AL_list},
        "results": results,
        "provenance": provenance
    }, f, indent=2)

pd.DataFrame(rows).to_csv(out_csv, index=False)

# Print concise human-readable summary
print("Fe–Al CALPHAD equilibrium at T = 973 K, P = 101325 Pa")
print(f"Database: {TDB_PATH}")
print(f"Phases considered: {len(relevant_phases)}")
print("")
for entry in results:
    x = entry['X_AL']
    phs = entry['phases']
    print(f"X_AL = {x:.2f}")
    if not phs:
        print("  (no phases above 1e-3 fraction)")
    for p in phs:
        print(f"  {p['phase']}: {p['mole_fraction']:.4f}")
    print("")

# Specific check for Al2Fe-like phases
any_al2fe_like = any(al2fe_like_present_by_x.values())
print("Al2Fe-like phase presence by composition (>=1e-3 fraction):")
for x in X_AL_list:
    yn = al2fe_like_present_by_x[float(x)]
    print(f"  X_AL = {x:.2f}: {'YES' if yn else 'NO'}")

# Also print a single-line JSON summary for quick parsing if needed
summary = {
    "al2fe_like_present": al2fe_like_present_by_x,
    "outputs": {"json": out_json, "csv": out_csv}
}
print(json.dumps(summary))
