#!/usr/bin/env python3
"""
CHGNet-based calculation of formation energies for nitrogen substitutional defects in GaAs:
 - N_Ga: N substitutes a Ga site
 - N_As: N substitutes an As site

Methodology (neutral defects, 0 K enthalpy approximation):
E_f(N_Ga) = E_tot(GaAs:N_Ga) - E_tot(GaAs) + mu_Ga - mu_N
E_f(N_As) = E_tot(GaAs:N_As) - E_tot(GaAs) + mu_As - mu_N

Chemical potentials chosen to be consistent with GaAs stability (As-rich) and N bounded by GaN stability:
 - mu_As = E(As_bulk)/atom
 - mu_Ga = mu_GaAs - mu_As,  where mu_GaAs = E(GaAs_bulk)/f.u.
 - mu_N  = mu_GaN - mu_Ga,   where mu_GaN  = E(GaN_bulk)/f.u.

Outputs:
 - results.json and summary.csv with energies and formation energies
 - final POSCAR-like VASP files for perfect and defect supercells and bulk references

Notes:
 - Loads CHGNet once; reuses calculator for all relaxations
 - Two-stage relaxation: positions-only then hydrostatic cell relax
 - STP temperature effects on chemical potentials are neglected (common for MLIP/DFT defect energetics)
"""
from __future__ import annotations

import os
import json
import math
import random
from dataclasses import dataclass
from typing import Dict, Tuple, List

import numpy as np

from dotenv import load_dotenv

from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ase import Atoms
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter

from chgnet.model import CHGNet, CHGNetCalculator


# ------------------------------ Utilities ------------------------------ #

def set_seeds(seed: int = 42) -> None:
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def ensure_env() -> None:
    load_dotenv()
    if not os.getenv("MP_API_KEY"):
        raise RuntimeError(
            "MP_API_KEY not found. Please set it in .env (MP_API_KEY=...) to query Materials Project."
        )


def get_best_structure(mpr: MPRester, formula: str) -> Structure:
    """
    Fetch the most stable structure for a given formula from Materials Project.
    Returns a pymatgen Structure.
    """
    docs = mpr.summary.search(
        formula=formula,
        fields=[
            "material_id",
            "formula_pretty",
            "energy_above_hull",
            "is_stable",
            "structure",
        ],
    )
    if not docs:
        raise RuntimeError(f"No MP entries found for {formula}")

    # Prefer stable, lowest energy-above-hull
    docs = sorted(docs, key=lambda d: (not d.is_stable, d.energy_above_hull))
    best = docs[0]
    if not hasattr(best, "structure") or best.structure is None:
        raise RuntimeError(f"MP entry for {formula} lacks a structure field")
    return best.structure


def to_conventional_supercell(struct: Structure, reps=(2, 2, 2)) -> Structure:
    """
    Convert to conventional standard structure (if available), then make a supercell.
    """
    try:
        sga = SpacegroupAnalyzer(struct, symprec=1e-3, angle_tolerance=5)
        conv = sga.get_conventional_standard_structure()
    except Exception:
        conv = struct.copy()
    sc = conv.copy()
    sc.make_supercell(reps)
    return sc


def write_poscar_from_atoms(atoms: Atoms, filename: str) -> None:
    ase_write(filename, atoms, format="vasp")


def atoms_num_formula_units(atoms: Atoms, formula: str) -> int:
    """
    For binary formula like 'GaAs' or 'GaN', return the number of formula units in the Atoms object.
    """
    from collections import Counter

    counts = Counter(atoms.get_chemical_symbols())
    # Parse formula composition (assume simple 1:1 like GaAs, GaN)
    # For generality, map symbol to required stoichiometry via pymatgen Composition
    from pymatgen.core import Composition

    comp = Composition(formula)
    # Normalize to integer coeffs
    elems = list(comp.as_dict().keys())
    ratios = [comp[el] for el in elems]
    # Determine integer multiplier to compare
    # For GaAs and GaN this is 1:1; robustly handle fractional by scaling
    lcm_denom = 1
    for r in ratios:
        frac = math.modf(r)[0]
        if abs(frac) > 1e-8:
            # scale denominator by 10 until integer-ish
            d = 1
            while abs(r * d - round(r * d)) > 1e-8 and d < 1000:
                d *= 10
            lcm_denom = max(lcm_denom, d)
    int_counts = [int(round(r * lcm_denom)) for r in ratios]
    # Compute formula units limited by minimum element ratio in the atoms object
    fus = []
    for el, need in zip(elems, int_counts):
        present = counts.get(el, 0)
        if need == 0:
            continue
        fus.append(present * lcm_denom // need)
    if not fus:
        raise ValueError(f"Could not determine formula units for {formula}")
    return int(min(fus))


@dataclass
class RelaxResult:
    name: str
    energy: float  # eV
    n_atoms: int
    atoms: Atoms


def check_finite(atoms: Atoms, label: str) -> None:
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not np.isfinite(e) or not np.isfinite(f).all():
        raise RuntimeError(f"Non-finite energy/forces for {label}")


def relax_structure(
    atoms: Atoms,
    calc: CHGNetCalculator,
    prefix: str,
    fmax_coarse: float = 0.2,
    fmax_fine: float = 0.05,
    max_steps: int = 300,
) -> RelaxResult:
    atoms = atoms.copy()
    atoms.set_calculator(calc)

    # Stage 1: positions-only relaxation
    log1 = f"{prefix}_relax1.log"
    opt1 = BFGS(atoms, logfile=log1)
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)
    check_finite(atoms, f"{prefix} after stage 1")

    # Stage 2: hydrostatic cell + positions relaxation
    log2 = f"{prefix}_relax2.log"
    try:
        ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
        opt2 = BFGS(ucf, logfile=log2)
        opt2.run(fmax=fmax_fine, steps=max_steps)
        # forces from filtered object are not directly accessible; re-evaluate
        check_finite(atoms, f"{prefix} after stage 2")
    except Exception as exc:
        # Fallback: positions-only fine relaxation
        with open(log2, "a") as fh:
            fh.write(f"\nUnitCellFilter failed with {exc}; falling back to positions-only fine relax.\n")
        opt2b = BFGS(atoms, logfile=log2)
        opt2b.run(fmax=fmax_fine, steps=max_steps)
        check_finite(atoms, f"{prefix} after stage 2 fallback")

    energy = atoms.get_potential_energy()
    return RelaxResult(prefix, float(energy), len(atoms), atoms)


def make_substitution(atoms: Atoms, from_symbol: str, to_symbol: str, prefix: str) -> Tuple[Atoms, int]:
    """Replace the first occurrence of from_symbol by to_symbol; return new Atoms and index replaced."""
    new_atoms = atoms.copy()
    symbols = new_atoms.get_chemical_symbols()
    try:
        idx = symbols.index(from_symbol)
    except ValueError:
        raise RuntimeError(f"No site with symbol {from_symbol} found in {prefix}")
    new_atoms[idx].symbol = to_symbol
    return new_atoms, idx


def main():
    set_seeds(42)
    ensure_env()

    # Load CHGNet once
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    adaptor = AseAtomsAdaptor()

    results: Dict[str, RelaxResult] = {}

    with MPRester(api_key=os.getenv("MP_API_KEY")) as mpr:
        # -------- Host GaAs perfect supercell --------
        gaas_struct = get_best_structure(mpr, "GaAs")
        gaas_sc = to_conventional_supercell(gaas_struct, reps=(2, 2, 2))
        gaas_atoms = adaptor.get_atoms(gaas_sc)
        write_poscar_from_atoms(gaas_atoms, "gaas_2x2x2_initial.vasp")
        res_gaias = relax_structure(gaas_atoms, calc, prefix="gaas_2x2x2")
        write_poscar_from_atoms(res_gaias.atoms, "gaas_2x2x2_relaxed.vasp")
        results[res_gaias.name] = res_gaias

        # Count formula units in the perfect supercell
        n_fu_gaas = atoms_num_formula_units(res_gaias.atoms, "GaAs")
        mu_gaas = res_gaias.energy / n_fu_gaas  # eV per GaAs formula unit

        # -------- Defects: N at As site (isovalent) --------
        nas_atoms, idx_as = make_substitution(res_gaias.atoms, from_symbol="As", to_symbol="N", prefix="N_As")
        write_poscar_from_atoms(nas_atoms, "gaas_2x2x2_N_As_initial.vasp")
        res_nas = relax_structure(nas_atoms, calc, prefix="gaas_2x2x2_N_As")
        write_poscar_from_atoms(res_nas.atoms, "gaas_2x2x2_N_As_relaxed.vasp")
        results[res_nas.name] = res_nas

        # -------- Defects: N at Ga site (heterovalent) --------
        nga_atoms, idx_ga = make_substitution(res_gaias.atoms, from_symbol="Ga", to_symbol="N", prefix="N_Ga")
        write_poscar_from_atoms(nga_atoms, "gaas_2x2x2_N_Ga_initial.vasp")
        res_nga = relax_structure(nga_atoms, calc, prefix="gaas_2x2x2_N_Ga")
        write_poscar_from_atoms(res_nga.atoms, "gaas_2x2x2_N_Ga_relaxed.vasp")
        results[res_nga.name] = res_nga

        # -------- Bulk references for chemical potentials --------
        # As bulk (grey As, A7)
        as_struct = get_best_structure(mpr, "As")
        as_atoms = adaptor.get_atoms(as_struct)
        res_as = relax_structure(as_atoms, calc, prefix="as_bulk")
        write_poscar_from_atoms(res_as.atoms, "as_bulk_relaxed.vasp")
        results[res_as.name] = res_as
        mu_as = res_as.energy / res_as.n_atoms  # eV/atom

        # Ga bulk (alpha-Ga)
        ga_struct = get_best_structure(mpr, "Ga")
        ga_atoms = adaptor.get_atoms(ga_struct)
        res_ga = relax_structure(ga_atoms, calc, prefix="ga_bulk")
        write_poscar_from_atoms(res_ga.atoms, "ga_bulk_relaxed.vasp")
        results[res_ga.name] = res_ga
        mu_ga_bulk = res_ga.energy / res_ga.n_atoms  # eV/atom

        # GaN bulk (wurtzite preferred)
        gan_struct = get_best_structure(mpr, "GaN")
        gan_atoms = adaptor.get_atoms(gan_struct)
        res_gan = relax_structure(gan_atoms, calc, prefix="gan_bulk")
        write_poscar_from_atoms(res_gan.atoms, "gan_bulk_relaxed.vasp")
        results[res_gan.name] = res_gan
        n_fu_gan = atoms_num_formula_units(res_gan.atoms, "GaN")
        mu_gan = res_gan.energy / n_fu_gan  # eV per GaN f.u.

        # As-rich condition for GaAs stability; define mu_Ga from GaAs
        mu_ga = mu_gaas - mu_as  # eV/atom
        # N chemical potential from GaN stability
        mu_n = mu_gan - mu_ga  # eV/atom of N

        # -------- Formation energies (neutral defects) --------
        e_perfect = res_gaias.energy
        e_n_as = res_nas.energy
        e_n_ga = res_nga.energy

        ef_n_as = (e_n_as - e_perfect) + (mu_as - mu_n)
        ef_n_ga = (e_n_ga - e_perfect) + (mu_ga - mu_n)
        delta = ef_n_ga - ef_n_as

        # Determine stability
        more_stable = "N_As" if ef_n_as < ef_n_ga else "N_Ga"

    # -------- Persist results --------
    summary = {
        "E_f_eV_N_Ga": ef_n_ga,
        "E_f_eV_N_As": ef_n_as,
        "Delta_E_f_eV": delta,
        "more_stable_site": more_stable,
        "notes": {
            "mu_As": mu_as,
            "mu_Ga": mu_ga,
            "mu_N": mu_n,
            "mu_GaAs": mu_gaas,
            "mu_GaN": mu_gan,
            "assumptions": "As-rich for GaAs; mu_N from GaN stability; neutral defects; 0 K enthalpy approximation",
        },
    }

    with open("results.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Write a small CSV summary
    with open("summary.csv", "w") as f:
        f.write("label,energy_eV,n_atoms\n")
        for k, v in results.items():
            f.write(f"{k},{v.energy},{v.n_atoms}\n")
        f.write("\n")
        f.write("quantity,value\n")
        f.write(f"E_f(N_Ga) [eV],{ef_n_ga}\n")
        f.write(f"E_f(N_As) [eV],{ef_n_as}\n")
        f.write(f"Delta E_f [eV],{delta}\n")
        f.write(f"more_stable,{more_stable}\n")

    # Human-readable printout
    print("Formation Energies (neutral defects, eV):")
    print(f"  E_f(N_Ga) = {ef_n_ga:.6f} eV")
    print(f"  E_f(N_As) = {ef_n_as:.6f} eV")
    print(f"  Delta     = {delta:.6f} eV  [E_f(N_Ga) - E_f(N_As)]")
    print(f"  More stable site: {more_stable}")


if __name__ == "__main__":
    main()

