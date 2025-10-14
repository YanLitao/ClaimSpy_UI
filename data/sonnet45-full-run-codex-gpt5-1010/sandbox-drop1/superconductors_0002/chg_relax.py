#!/usr/bin/env python3
"""
CHGNet-based stability comparison for Zr3Sb A15 vs Ni3P at ~0 GPa.

This script:
  1) Loads MP_API_KEY from .env
  2) Fetches Zr3Sb structures (A15 and Ni3P types) from Materials Project
  3) Fetches elemental references (Zr, Sb) and stable Zr–Sb phases (for hull)
  4) Runs CHGNet relaxations (two-stage) for all relevant structures
  5) Computes CHGNet formation energies, builds a binary convex hull, and
     reports above-hull energies at 25 at.% Sb (x_Sb=0.25)
  6) Writes results.json, summary.csv, and POSCARs for initial/final structures.

Notes:
  - Uses a single CHGNet instance and reuses the calculator across relaxations.
  - Saves logs but never blocks for input; no GUI.
  - Prints a concise human-readable summary.
"""

from __future__ import annotations

import json
import math
import os
import sys
import csv
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from dotenv import load_dotenv

# MP API
from mp_api.client import MPRester

# Pymatgen/ASE/CHGNet
from pymatgen.core import Structure
from pymatgen.core.composition import Composition
from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

from ase import Atoms
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator


# ----------------------------- Utilities -----------------------------------


def set_seed(seed: int = 42) -> None:
    import random
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass


def write_poscar(struct: Structure, path: str) -> None:
    Poscar(struct).write_file(path)


def composition_tuple(struct: Structure) -> Tuple[int, int]:
    comp = struct.composition.get_el_amt_dict()
    zr = int(round(comp.get("Zr", 0)))
    sb = int(round(comp.get("Sb", 0)))
    if zr + sb == 0:
        # Non Zr/Sb - not expected here, but guard.
        zr = int(round(struct.composition["Zr"])) if "Zr" in struct.composition else 0
        sb = int(round(struct.composition["Sb"])) if "Sb" in struct.composition else 0
    return zr, sb


def sb_fraction(struct: Structure) -> float:
    zr, sb = composition_tuple(struct)
    total = zr + sb
    return sb / total if total > 0 else 0.0


def label_from_doc(doc) -> str:
    # Prefer formula_pretty; fall back to composition_reduced or material_id
    if getattr(doc, "formula_pretty", None):
        return doc.formula_pretty
    if getattr(doc, "composition_reduced", None):
        try:
            return Composition(doc.composition_reduced).reduced_formula
        except Exception:
            pass
    return getattr(doc, "material_id", "unknown")


def structure_type_hint(doc) -> Optional[str]:
    # Use spacegroup symbol
    sym = getattr(doc, "symmetry", None)
    if sym and getattr(sym, "symbol", None):
        sg = sym.symbol
        if sg == "Pm-3n":
            return "A15"
        if sg == "I-4":
            return "NI3P"
    return None


def safe_min_distance_scale(atoms: Atoms, min_allow: float = 0.7) -> None:
    # Very conservative fix: if any interatomic distance < min_allow Å, scale cell isotropically by 5% until safe.
    # This should be rare for MP-derived structures.
    for _ in range(5):
        try:
            dmin = atoms.get_distance(0, 1, mic=True) if len(atoms) > 1 else 1e3
            # That's only one pair; do a more thorough check when small system
            if len(atoms) <= 64:
                dmin = 1e9
                pos = atoms.get_positions()
                for i in range(len(atoms)):
                    for j in range(i + 1, len(atoms)):
                        d = atoms.get_distance(i, j, mic=True)
                        if d < dmin:
                            dmin = d
            if dmin >= min_allow:
                return
            # Scale
            cell = atoms.get_cell()
            atoms.set_cell(cell * 1.05, scale_atoms=True)
        except Exception:
            break


@dataclass
class RelaxResult:
    label: str
    mp_id: Optional[str]
    structure_type: Optional[str]
    spacegroup: Optional[str]
    x_sb: float
    natoms: int
    energy_per_atom: float
    formation_energy_per_atom: Optional[float] = None
    above_hull_eV_per_atom: Optional[float] = None
    init_poscar: Optional[str] = None
    final_poscar: Optional[str] = None


def relax_with_chgnet(
    atoms: Atoms,
    calc: CHGNetCalculator,
    label: str,
    workdir: str,
    fmax_coarse: float = 0.2,
    fmax_fine: float = 0.05,
    max_steps: int = 300,
) -> Tuple[Atoms, float]:
    # Assign calculator once
    atoms.set_calculator(calc)

    # Stage 1: positions-only
    opt1 = BFGS(atoms, logfile=os.path.join(workdir, f"{label}_relax1.log"))
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)

    # Stage 2: limited cell relax (hydrostatic strain)
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile=os.path.join(workdir, f"{label}_relax2.log"))
    opt2.run(fmax=fmax_fine, steps=max_steps)

    epa = atoms.get_potential_energy() / len(atoms)
    return atoms, float(epa)


def build_lower_hull(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Monotone chain lower convex hull for (x, y) with x in [0,1]."""
    pts = sorted(points, key=lambda t: (t[0], t[1]))
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    return lower


def hull_energy_at(x: float, lower_hull: List[Tuple[float, float]]) -> Tuple[float, Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Linear interpolation of hull energy at composition x (0..1). Returns energy and the bracketing segment."""
    if x <= lower_hull[0][0]:
        return lower_hull[0][1], (lower_hull[0], lower_hull[0])
    if x >= lower_hull[-1][0]:
        return lower_hull[-1][1], (lower_hull[-1], lower_hull[-1])
    for i in range(1, len(lower_hull)):
        x0, y0 = lower_hull[i - 1]
        x1, y1 = lower_hull[i]
        if x0 <= x <= x1:
            if abs(x1 - x0) < 1e-12:
                return min(y0, y1), ((x0, y0), (x1, y1))
            t = (x - x0) / (x1 - x0)
            y = y0 * (1 - t) + y1 * t
            return y, ((x0, y0), (x1, y1))
    # Fallback (should not happen)
    return lower_hull[-1][1], (lower_hull[-2], lower_hull[-1])


def main() -> None:
    set_seed(42)

    workdir = os.getcwd()
    load_dotenv(dotenv_path=os.path.join(workdir, ".env"), override=False)
    api_key = os.environ.get("MP_API_KEY")
    if not api_key:
        print("ERROR: MP_API_KEY not found in .env", file=sys.stderr)
        sys.exit(1)

    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)
    adaptor = AseAtomsAdaptor()

    zr3sb_docs: List = []
    elemental_docs: Dict[str, Optional[object]] = {"Zr": None, "Sb": None}
    zr_sb_stable_docs: List = []

    with MPRester(api_key) as mpr:
        # Fetch Zr3Sb candidates (all structures); filter locally by type
        zr3sb_docs = mpr.summary.search(
            formula=["Zr3Sb"],
            fields=[
                "material_id",
                "formula_pretty",
                "energy_per_atom",
                "structure",
                "symmetry",
                "deprecated",
                "energy_above_hull",
            ],
        )

        # Elemental references
        for el in ["Zr", "Sb"]:
            docs = mpr.summary.search(
                formula=[el],
                is_stable=True,
                fields=[
                    "material_id",
                    "formula_pretty",
                    "structure",
                    "symmetry",
                ],
            )
            elemental_docs[el] = docs[0] if docs else None

        # Stable Zr–Sb phases to approximate hull
        zr_sb_stable_docs = mpr.summary.search(
            chemsys=["Zr-Sb"],
            num_elements=2,
            is_stable=True,
            fields=[
                "material_id",
                "formula_pretty",
                "structure",
                "symmetry",
                "energy_above_hull",
            ],
        )

    # Validate elemental docs
    for el, doc in elemental_docs.items():
        if doc is None:
            print(f"ERROR: Could not fetch elemental reference for {el}", file=sys.stderr)
            sys.exit(2)

    # Select Zr3Sb A15 and Ni3P candidates
    a15_doc = None
    ni3p_doc = None
    for doc in zr3sb_docs:
        hint = structure_type_hint(doc)
        if hint == "A15" and a15_doc is None:
            a15_doc = doc
        if hint == "NI3P" and ni3p_doc is None:
            ni3p_doc = doc

    # If any missing, try fallback: check by spacegroup directly among docs
    if a15_doc is None:
        for doc in zr3sb_docs:
            sym = getattr(doc, "symmetry", None)
            if sym and getattr(sym, "symbol", None) == "Pm-3n":
                a15_doc = doc
                break
    if ni3p_doc is None:
        for doc in zr3sb_docs:
            sym = getattr(doc, "symmetry", None)
            if sym and getattr(sym, "symbol", None) == "I-4":
                ni3p_doc = doc
                break

    # If still missing, try broader search in Zr–Sb system for A15/Ni3P with formula Zr3Sb
    # (should be rare; keeping code compact)
    if (a15_doc is None) or (ni3p_doc is None):
        with MPRester(api_key) as mpr:
            docs = mpr.summary.search(
                chemsys=["Zr-Sb"],
                fields=[
                    "material_id",
                    "formula_pretty",
                    "structure",
                    "symmetry",
                ],
            )
        for doc in docs:
            # Filter exact composition Zr3Sb
            if label_from_doc(doc) != "Zr3Sb":
                continue
            hint = structure_type_hint(doc)
            if hint == "A15" and a15_doc is None:
                a15_doc = doc
            if hint == "NI3P" and ni3p_doc is None:
                ni3p_doc = doc

    if a15_doc is None or ni3p_doc is None:
        # Attempt to construct A15 by species substitution from a canonical A15 (e.g., Cr3Si)
        if a15_doc is None:
            with MPRester(api_key) as mpr:
                cands = mpr.summary.search(
                    formula=["Cr3Si", "Nb3Sn", "V3Si", "Mo3Si"],
                    fields=["material_id", "formula_pretty", "structure", "symmetry"],
                )
            # Pick first with Pm-3n
            proto_doc = None
            for d in cands:
                sym = getattr(d, "symmetry", None)
                if sym and getattr(sym, "symbol", None) == "Pm-3n":
                    proto_doc = d
                    break
            if proto_doc is None:
                print(
                    "ERROR: Could not locate any canonical A15 prototype to substitute from.",
                    file=sys.stderr,
                )
                sys.exit(3)
            # Substitute A->Zr, B->Sb
            prot = proto_doc.structure.copy()
            species = set([str(sp) for sp in prot.composition.elements])
            # Heuristic mapping: species sorted alphabetically -> (A,B)
            # Ensure stoichiometry A3B
            comp = prot.composition.get_el_amt_dict()
            # Identify majority species as A (3), minority as B (1)
            items = sorted(comp.items(), key=lambda kv: kv[1], reverse=True)
            if len(items) < 2:
                print("ERROR: Prototype does not have two species.", file=sys.stderr)
                sys.exit(3)
            A_sp = items[0][0]
            B_sp = items[1][0]
            sub_map = {A_sp: "Zr", B_sp: "Sb"}
            try:
                prot_sub = prot.copy()
                prot_sub.replace_species(sub_map)
            except Exception as exc:
                print(f"ERROR: Substitution failed: {exc}", file=sys.stderr)
                sys.exit(3)

            class Dummy:
                pass

            a15_doc = Dummy()
            a15_doc.material_id = None
            a15_doc.formula_pretty = "Zr3Sb"
            a15_doc.structure = prot_sub
            class Sym:
                def __init__(self):
                    self.symbol = "Pm-3n"
            a15_doc.symmetry = Sym()

        # If Ni3P still missing, fail (rarer)
        if ni3p_doc is None:
            print("ERROR: Could not locate Zr3Sb Ni3P prototype from MP.", file=sys.stderr)
            sys.exit(3)

    # Prepare all structures to relax
    targets: List[Tuple[str, Optional[str], Optional[str], Optional[str], Structure]] = []

    def add_target(doc, tag: str):
        stype = structure_type_hint(doc)
        sg = getattr(getattr(doc, "symmetry", None), "symbol", None)
        targets.append((tag, getattr(doc, "material_id", None), stype, sg, doc.structure))

    add_target(a15_doc, "Zr3Sb_A15")
    add_target(ni3p_doc, "Zr3Sb_Ni3P")
    # elementals
    for el in ["Zr", "Sb"]:
        doc = elemental_docs[el]
        sg = getattr(getattr(doc, "symmetry", None), "symbol", None)
        targets.append((el, getattr(doc, "material_id", None), el, sg, doc.structure))
    # Zr–Sb stable phases (for hull)
    # Deduplicate by reduced formula to avoid multiple polymorphs
    seen_rf = set()
    for doc in zr_sb_stable_docs:
        rf = label_from_doc(doc)
        if rf in seen_rf:
            continue
        seen_rf.add(rf)
        # Skip Zr3Sb A15/Ni3P here; those are already added explicitly
        if rf == "Zr3Sb":
            continue
        stype = structure_type_hint(doc)
        sg = getattr(getattr(doc, "symmetry", None), "symbol", None)
        targets.append((rf, getattr(doc, "material_id", None), stype, sg, doc.structure))

    # Relax all targets
    results: List[RelaxResult] = []
    e_ref: Dict[str, float] = {}

    for tag, mp_id, stype, sg, pmg_struct in targets:
        label = tag
        # Save initial POSCARs for key ones
        init_poscar_path = os.path.join(workdir, f"POSCAR_{label}.vasp")
        try:
            write_poscar(pmg_struct, init_poscar_path)
        except Exception:
            init_poscar_path = None

        atoms = adaptor.get_atoms(pmg_struct)
        safe_min_distance_scale(atoms)

        atoms_relaxed, epa = relax_with_chgnet(atoms, calc, label, workdir)

        # Save final structure back to POSCAR
        try:
            final_struct: Structure = adaptor.get_structure(atoms_relaxed)
            final_poscar_path = os.path.join(workdir, f"POSCAR_{label}_relaxed.vasp")
            write_poscar(final_struct, final_poscar_path)
        except Exception:
            final_poscar_path = None

        xsb = sb_fraction(pmg_struct)
        res = RelaxResult(
            label=label,
            mp_id=mp_id,
            structure_type=stype,
            spacegroup=sg,
            x_sb=xsb,
            natoms=len(atoms_relaxed),
            energy_per_atom=epa,
            init_poscar=init_poscar_path,
            final_poscar=final_poscar_path,
        )

        results.append(res)

        if label in ("Zr", "Sb"):
            e_ref[label] = epa

    if not ("Zr" in e_ref and "Sb" in e_ref):
        print("ERROR: Missing elemental reference energies after relaxation.", file=sys.stderr)
        sys.exit(4)

    # Compute formation energies
    def formation_energy_per_atom(res: RelaxResult) -> float:
        if res.label in ("Zr", "Sb"):
            return 0.0
        # Use reduced composition from label if possible; fall back to x_sb
        # For binary Zr–Sb, we can derive stoichiometry from initial structure
        # (using the init POSCAR path is not reliable). Instead, use x_sb and total atoms in reduced formula.
        # We'll take composition from the final structure positions by mapping back to initial pmg struct if needed.
        # However, we retained the original pmg_struct composition via x_sb; also need total atom counts.
        # Safer: parse from label if in form Zr[a]Sb[b]. Otherwise, infer from x_sb by rounding to small integers.
        zr = None
        sb = None
        # Try parse label
        if res.label.startswith("Zr") and "Sb" in res.label:
            try:
                tail = res.label[len("Zr"):]
                zs = ""
                i = 0
                while i < len(tail) and tail[i].isdigit():
                    zs += tail[i]
                    i += 1
                zr = int(zs) if zs else 1
                assert tail[i:i+2] == "Sb"
                i += 2
                ss = ""
                while i < len(tail) and tail[i].isdigit():
                    ss += tail[i]
                    i += 1
                sb = int(ss) if ss else 1
            except Exception:
                zr = None
                sb = None
        if zr is None or sb is None:
            # Infer small integers up to 12 atoms
            x = res.x_sb
            best = (None, None, 1e9)
            for tot in range(2, 13):
                s = round(x * tot)
                z = tot - s
                if s <= 0 or z <= 0:
                    continue
                if abs(s / tot - x) < best[2]:
                    best = (z, s, abs(s / tot - x))
            zr, sb, _ = best
        m = int(zr)
        n = int(sb)
        e_tot = res.energy_per_atom * (m + n)
        e_form = (e_tot - m * e_ref["Zr"] - n * e_ref["Sb"]) / (m + n)
        return float(e_form)

    for r in results:
        r.formation_energy_per_atom = formation_energy_per_atom(r)

    # Build binary hull points: include elements at (x=0,0) and (x=1,0)
    points: List[Tuple[float, float]] = []
    # Consolidate by composition (x_sb) keep the lowest formation energy at each x to avoid duplicates
    best_at_x: Dict[float, float] = {}
    for r in results:
        if r.label in ("Zr", "Sb"):
            continue
        x = r.x_sb
        y = r.formation_energy_per_atom
        if math.isfinite(x) and math.isfinite(y):
            if (x not in best_at_x) or (y < best_at_x[x]):
                best_at_x[x] = y
    # Add elements
    points.append((0.0, 0.0))
    points.append((1.0, 0.0))
    for x, y in sorted(best_at_x.items()):
        points.append((x, y))

    lower = build_lower_hull(points)
    x_target = 0.25
    hull_e, seg = hull_energy_at(x_target, lower)

    # Assign above-hull for all
    for r in results:
        if r.label in ("Zr", "Sb"):
            r.above_hull_eV_per_atom = 0.0
            continue
        x = r.x_sb
        he, _ = hull_energy_at(x, lower)
        r.above_hull_eV_per_atom = r.formation_energy_per_atom - he

    # Extract A15 and Ni3P results
    a15_res = next((r for r in results if r.label == "Zr3Sb_A15"), None)
    ni3p_res = next((r for r in results if r.label == "Zr3Sb_Ni3P"), None)
    if a15_res is None or ni3p_res is None:
        print("ERROR: Missing prototype results.", file=sys.stderr)
        sys.exit(5)

    # Decide ground state among prototypes
    ground_label = "A15" if a15_res.energy_per_atom < ni3p_res.energy_per_atom else "Ni3P"
    delta_e = abs(a15_res.energy_per_atom - ni3p_res.energy_per_atom)

    # Identify competing phases (tie-line endpoints) near x=0.25
    # Map hull segment endpoints back to nearest computed phases (by x,y)
    def find_nearest_phase(xy: Tuple[float, float]) -> Optional[RelaxResult]:
        x, y = xy
        best = None
        bestd = 1e9
        for r in results:
            if r.label in ("Zr", "Sb"):
                continue
            d = abs(r.x_sb - x) + abs(r.formation_energy_per_atom - y)
            if d < bestd:
                bestd = d
                best = r
        return best

    left, right = seg
    left_phase = find_nearest_phase(left)
    right_phase = find_nearest_phase(right)

    # Prepare outputs
    results_json = {
        "references": {
            "Zr_energy_per_atom": e_ref["Zr"],
            "Sb_energy_per_atom": e_ref["Sb"],
        },
        "prototypes": {
            "A15": {
                "energy_per_atom": a15_res.energy_per_atom,
                "formation_energy_per_atom": a15_res.formation_energy_per_atom,
                "above_hull_eV_per_atom": a15_res.above_hull_eV_per_atom,
                "mp_id": a15_res.mp_id,
                "spacegroup": a15_res.spacegroup,
            },
            "Ni3P": {
                "energy_per_atom": ni3p_res.energy_per_atom,
                "formation_energy_per_atom": ni3p_res.formation_energy_per_atom,
                "above_hull_eV_per_atom": ni3p_res.above_hull_eV_per_atom,
                "mp_id": ni3p_res.mp_id,
                "spacegroup": ni3p_res.spacegroup,
            },
        },
        "ground_state_among_prototypes": ground_label,
        "energy_difference_eV_per_atom": delta_e,
        "hull_at_xSb_0.25": {
            "hull_energy_eV_per_atom": hull_e,
            "segment": {"left": left, "right": right},
            "competing_left": {
                "label": getattr(left_phase, "label", None),
                "mp_id": getattr(left_phase, "mp_id", None),
                "x_sb": getattr(left_phase, "x_sb", None),
                "E_form": getattr(left_phase, "formation_energy_per_atom", None),
            },
            "competing_right": {
                "label": getattr(right_phase, "label", None),
                "mp_id": getattr(right_phase, "mp_id", None),
                "x_sb": getattr(right_phase, "x_sb", None),
                "E_form": getattr(right_phase, "formation_energy_per_atom", None),
            },
        },
        "all_results": [
            {
                "label": r.label,
                "mp_id": r.mp_id,
                "structure_type": r.structure_type,
                "spacegroup": r.spacegroup,
                "x_sb": r.x_sb,
                "natoms": r.natoms,
                "energy_per_atom": r.energy_per_atom,
                "formation_energy_per_atom": r.formation_energy_per_atom,
                "above_hull_eV_per_atom": r.above_hull_eV_per_atom,
                "init_poscar": r.init_poscar,
                "final_poscar": r.final_poscar,
            }
            for r in results
        ],
    }

    json_path = os.path.join(workdir, "results.json")
    with open(json_path, "w") as f:
        json.dump(results_json, f, indent=2)

    # CSV summary
    csv_path = os.path.join(workdir, "summary.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "label",
                "mp_id",
                "structure_type",
                "spacegroup",
                "x_sb",
                "energy_per_atom(eV)",
                "formation_energy_per_atom(eV)",
                "above_hull_eV_per_atom(eV)",
            ]
        )
        for r in results:
            w.writerow(
                [
                    r.label,
                    r.mp_id,
                    r.structure_type,
                    r.spacegroup,
                    f"{r.x_sb:.6f}",
                    f"{r.energy_per_atom:.6f}",
                    f"{r.formation_energy_per_atom:.6f}",
                    f"{r.above_hull_eV_per_atom:.6f}" if r.above_hull_eV_per_atom is not None else "",
                ]
            )

    # Human summary printout
    print("CHGNet Zr3Sb prototype comparison at ~0 GPa (static relax):")
    print(f"- A15:    E/atom = {a15_res.energy_per_atom:.6f} eV, E_form = {a15_res.formation_energy_per_atom:.6f} eV/atom, E_above_hull = {a15_res.above_hull_eV_per_atom:.6f} eV/atom")
    print(f"- Ni3P:   E/atom = {ni3p_res.energy_per_atom:.6f} eV, E_form = {ni3p_res.formation_energy_per_atom:.6f} eV/atom, E_above_hull = {ni3p_res.above_hull_eV_per_atom:.6f} eV/atom")
    print(f"- Ground state among these: {ground_label} (ΔE = {delta_e:.6f} eV/atom)")
    is_a15_near = (a15_res.above_hull_eV_per_atom is not None) and (a15_res.above_hull_eV_per_atom <= 0.025)
    print(f"- A15 near hull (≤25 meV/atom)? {'YES' if is_a15_near else 'NO'} (E_above_hull={a15_res.above_hull_eV_per_atom:.6f} eV/atom)")
    ltxt = f"{getattr(left_phase, 'label', 'unknown')} (x_Sb={getattr(left_phase, 'x_sb', None)})"
    rtxt = f"{getattr(right_phase, 'label', 'unknown')} (x_Sb={getattr(right_phase, 'x_sb', None)})"
    print(f"- Hull at 25% Sb set by tie-line: {ltxt} ↔ {rtxt}; E_hull(0.25) = {hull_e:.6f} eV/atom")
    print(f"Outputs: results.json, summary.csv, POSCAR_* for initial/final structures")


if __name__ == "__main__":
    main()
