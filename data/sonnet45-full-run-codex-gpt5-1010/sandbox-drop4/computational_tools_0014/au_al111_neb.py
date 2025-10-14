#!/usr/bin/env python3
"""
Au adatom diffusion on Al(111) using CHGNet + ASE NEB

Workflow:
- Build Al(111) slab (pymatgen SlabGenerator; 3x3 in-plane, 4 layers, ~12 Å vacuum)
- Place Au adatom at two adjacent fcc hollow sites (AdsorbateSiteFinder)
- Convert to ASE; fix bottom two Al layers
- Relax endpoints (two-stage BFGS; positions-only)
- Run CI-NEB with IDPP interpolation
- Report forward barrier E_d = E_max - E_initial (eV)
- Save results: JSON, CSV, POSCARs, and NEB trajectory

Notes:
- CHGNet model is loaded exactly once and reused across calculators
- No interactive plotting; logs written to files
"""

import json
import os
import math
import numpy as np
from dataclasses import asdict, dataclass

# Ensure .env variables (e.g., MP_API_KEY) are loaded if present
try:
    from dotenv import load_dotenv
    load_dotenv()
except Exception:
    pass

from chgnet.model import CHGNet, CHGNetCalculator

from ase.build import fcc111, add_adsorbate
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

from ase.atoms import Atoms
from ase.io import write, Trajectory
from ase.constraints import FixAtoms
from ase.optimize import BFGS, FIRE
from ase.neb import NEB, NEBTools


np.random.seed(42)


@dataclass
class BarrierResult:
    barrier_eV: float
    barrier_method: str
    initial_energy_eV: float
    final_energy_eV: float
    max_energy_eV: float
    max_image_index: int
    n_images: int
    n_atoms: int
    slab_shape: str
    notes: str


def build_al111_slab(a0=4.05, nlayers=4, vacuum=12.0, size=(3, 3, 4)):
    # Build Al(111) slab using ASE with special adsorption site keywords supported
    slab = fcc111('Al', size=size, a=a0, vacuum=vacuum, orthogonal=False)
    return slab


def classify_hollow_sites_as_fcc_hcp(pmg_slab, hollow_sites, z_tol=0.4):
    # Cluster Al atom z-levels
    al_z = [site.coords[2] for site in pmg_slab.sites if site.species_string == 'Al']
    al_z_sorted = sorted(al_z)
    levels = []
    for z in al_z_sorted:
        if not levels or abs(z - levels[-1]) > z_tol:
            levels.append(z)
    if len(levels) < 3:
        # Not enough layers recognized to distinguish fcc vs hcp
        return {'fcc': [], 'hcp': []}
    z_top = max(al_z_sorted)
    # Build mapping from atom z to layer idx
    def layer_index(z):
        return int(np.argmin([abs(z - L) for L in levels]))
    # Determine top layer index
    top_index = layer_index(z_top)
    # Classify each hollow: find nearest underlying Al atom and check its layer
    from math import hypot
    out = {'fcc': [], 'hcp': []}
    al_positions = np.array([site.coords for site in pmg_slab.sites if site.species_string == 'Al'])
    al_layers = [layer_index(z) for z in al_positions[:, 2]]
    for pt in hollow_sites:
        # consider only atoms below the site
        below = al_positions[:, 2] < (pt[2] - 0.1)
        if not np.any(below):
            continue
        xy = np.array([pt[0], pt[1]])
        dels = al_positions[below][:, :2] - xy
        d2 = (dels ** 2).sum(axis=1)
        j = int(np.argmin(d2))
        idx_global = np.arange(len(al_positions))[below][j]
        lidx = al_layers[idx_global]
        # Relative layer index from the top: 0=top, -1=2nd-from-top etc.
        # Since we clustered ascending, compare to top_index
        rel = lidx - top_index
        # Top_index is the last in levels; but layer_index returns index in ascending order
        # So top_index should be len(levels)-1
        # Identify whether nearest atom is in the second layer (hcp) or third layer (fcc)
        if lidx == (len(levels) - 2):
            out['hcp'].append(pt)
        elif lidx == (len(levels) - 3):
            out['fcc'].append(pt)
    return out


def find_two_adjacent_fcc_sites_xy(slab: Atoms):
    adaptor = AseAtomsAdaptor()
    pmg_slab = adaptor.get_structure(slab)
    asf = AdsorbateSiteFinder(pmg_slab)
    sites = asf.find_adsorption_sites(symm_reduce=0)
    hollows = sites.get('hollow', [])
    classified = classify_hollow_sites_as_fcc_hcp(pmg_slab, hollows)
    fcc_sites = classified.get('fcc', [])
    if len(fcc_sites) < 2:
        # Fall back: pick two nearest hollows (may mix fcc/hcp)
        fcc_sites = hollows
    # Select top-surface fcc by highest z
    zmax = max(pt[2] for pt in fcc_sites)
    top = [pt for pt in fcc_sites if (zmax - pt[2]) < 1.0]
    if len(top) < 2:
        raise RuntimeError('Could not identify two top-surface fcc sites')
    # Choose central one and its nearest neighbor (distinct in-cell site)
    # Use 2D distances
    xy = np.array([[p[0], p[1]] for p in top])
    center = np.mean(slab.cell[:2], axis=0)
    i0 = int(np.argmin(((xy - center[:2]) ** 2).sum(axis=1)))
    d2 = ((xy - xy[i0]) ** 2).sum(axis=1)
    d2[i0] = 1e9
    i1 = int(np.argmin(d2))
    s1 = np.array([top[i0][0], top[i0][1]])
    s2 = np.array([top[i1][0], top[i1][1]])
    return s1, s2


def make_initial_final_with_au(slab: Atoms, height=2.2, base_offset=(1, 1)):
    # Create initial and final states with Au at two distinct in-cell fcc sites
    s1xy, s2xy = find_two_adjacent_fcc_sites_xy(slab)
    initial = slab.copy()
    final = slab.copy()
    add_adsorbate(initial, 'Au', height=height, position=(float(s1xy[0]), float(s1xy[1])))
    add_adsorbate(final, 'Au', height=height, position=(float(s2xy[0]), float(s2xy[1])))
    return initial, final


    # (pymatgen conversion helpers removed; using ASE builders directly)


def fix_bottom_layers(atoms: Atoms, n_layers_to_fix=2, z_tol=0.4):
    # Identify Al layers along z and fix the bottom n_layers_to_fix layers
    # Work only on Al atoms; Au remains free
    z_al = []
    for i, s in enumerate(atoms.get_chemical_symbols()):
        if s == "Al":
            z_al.append((i, atoms.positions[i, 2]))
    if not z_al:
        return
    # Sort by z
    z_al.sort(key=lambda x: x[1])
    levels = []
    for _, z in z_al:
        if not levels or abs(z - levels[-1]) > z_tol:
            levels.append(z)
    if len(levels) < n_layers_to_fix:
        # Fall back: fix half of the slab by z
        z_cut = np.median([z for _, z in z_al])
        mask = [(s == "Al" and atoms.positions[i, 2] < z_cut) for i, s in enumerate(atoms.get_chemical_symbols())]
    else:
        # Determine z boundary for the bottom n layers
        z_boundary = levels[n_layers_to_fix - 1] + z_tol / 2.0
        mask = [(s == "Al" and atoms.positions[i, 2] <= z_boundary) for i, s in enumerate(atoms.get_chemical_symbols())]
    atoms.set_constraint(FixAtoms(mask=mask))


def make_calc(chgnet_model):
    return CHGNetCalculator(chgnet_model)


def precheck_finite(atoms: Atoms):
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not (np.isfinite(e) and np.all(np.isfinite(f))):
        raise RuntimeError("Non-finite energy/forces encountered.")


def relax_positions(atoms: Atoms, calc: CHGNetCalculator, fmax_coarse=0.2, fmax_fine=0.05,
                    max_steps=300, logprefix="relax"):
    atoms.set_calculator(calc)
    precheck_finite(atoms)
    opt1 = BFGS(atoms, logfile=f"{logprefix}_1.log")
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)
    precheck_finite(atoms)
    opt2 = BFGS(atoms, logfile=f"{logprefix}_2.log")
    opt2.run(fmax=fmax_fine, steps=max_steps)
    precheck_finite(atoms)
    e = atoms.get_potential_energy()
    return e


def run_neb(initial: Atoms, final: Atoms, chgnet_model, n_images=5, fmax=0.05, max_steps=600,
            logprefix="neb"):
    # Prepare images
    images = [initial]
    for _ in range(n_images - 2):
        images.append(initial.copy())
    images.append(final)

    # Ensure constraints are applied to all images (clone from initial may already have them)
    for img in images:
        # Reapply fixation mask to be safe
        fix_bottom_layers(img, n_layers_to_fix=2)
        img.set_calculator(make_calc(chgnet_model))

    neb = NEB(images, climb=True)
    # Interpolate between initial and final; prefer IDPP if available
    try:
        neb.interpolate(method="idpp", mic=True)
    except TypeError:
        try:
            from ase.neb import idpp_interpolate
            idpp_interpolate(images)
        except Exception:
            neb.interpolate(mic=True)

    opt = BFGS(neb, logfile=f"{logprefix}.log")
    opt.run(fmax=fmax, steps=max_steps)

    # Collect energies after NEB
    energies = [img.get_potential_energy() for img in images]
    e0 = energies[0]
    emax = max(energies)
    imax = int(np.argmax(energies))

    # Also compute spline-refined barrier via NEBTools for comparison
    try:
        tools = NEBTools(images)
        # Returns (E_forward, E_backward), with reference to initial and final
        # Newer ASE: get_barrier() -> (Ef, Eb), accordingly
        Ef, Eb = tools.get_barrier()
        barrier = float(Ef)
        method = "CI-NEB (spline-refined)"
    except Exception:
        barrier = float(emax - e0)
        method = "CI-NEB (max-image)"

    return barrier, energies, imax, images


def main():
    # File naming
    results_json = "results.json"
    results_csv = "summary.csv"
    poscar_init = "POSCAR_initial.vasp"
    poscar_final = "POSCAR_final.vasp"
    poscar_ts = "POSCAR_ts_guess.vasp"
    neb_traj = "neb.traj"

    # 1) Build slab
    slab = build_al111_slab(a0=4.05, nlayers=4, vacuum=12.0, size=(3, 3, 4))

    # 2) Build initial and final with Au at adjacent fcc sites
    initial, final = make_initial_final_with_au(slab, height=2.2, base_offset=(1, 1))
    # Debug info: print Au positions (xy) before relaxation
    try:
        idx_init = [i for i, s in enumerate(initial.get_chemical_symbols()) if s == 'Au'][0]
        idx_final = [i for i, s in enumerate(final.get_chemical_symbols()) if s == 'Au'][0]
        pos_i = initial.positions[idx_init]
        pos_f = final.positions[idx_final]
        print(f"Initial Au xy: {pos_i[:2]}")
        print(f"Final   Au xy: {pos_f[:2]}")
    except Exception:
        pass
    for at in (initial, final):
        at.set_pbc([True, True, True])
        fix_bottom_layers(at, n_layers_to_fix=2)

    # 3) Load CHGNet once
    chgnet = CHGNet.load()

    # 4) Relax endpoints
    e_init = relax_positions(initial, make_calc(chgnet), logprefix="initial_relax")
    e_final = relax_positions(final, make_calc(chgnet), logprefix="final_relax")

    # Save endpoints
    write(poscar_init, initial, format="vasp")
    write(poscar_final, final, format="vasp")

    # 5) NEB: interpolate + optimize
    barrier, energies, imax, images = run_neb(initial, final, chgnet, n_images=7, fmax=0.05, max_steps=600, logprefix="neb")

    # Save a TS-guess POSCAR (max-energy image) and NEB trajectory
    try:
        write(poscar_ts, images[imax], format="vasp")
    except Exception:
        pass
    try:
        traj = Trajectory(neb_traj, "w")
        for img in images:
            traj.write(img)
        traj.close()
    except Exception:
        pass

    # 6) Summary and persistence
    n_atoms = len(initial)
    slab_shape = f"Al(111) 3x3x{4} + Au adatom"
    res = BarrierResult(
        barrier_eV=float(barrier),
        barrier_method="CI-NEB (IDPP init, BFGS, climb)",
        initial_energy_eV=float(e_init),
        final_energy_eV=float(e_final),
        max_energy_eV=float(max(energies)),
        max_image_index=int(imax),
        n_images=7,
        n_atoms=n_atoms,
        slab_shape=slab_shape,
        notes="Bottom two Al layers fixed; endpoints relaxed to fmax=0.05 eV/Å; NEB fmax=0.05 eV/Å.",
    )

    # Write JSON
    out = asdict(res)
    out["image_energies_eV"] = [float(x) for x in energies]
    with open(results_json, "w") as f:
        json.dump(out, f, indent=2)

    # Write CSV summary
    with open(results_csv, "w") as f:
        f.write("barrier_eV,initial_eV,final_eV,max_eV,max_image_index,n_images,n_atoms\n")
        f.write(
            f"{res.barrier_eV:.6f},{res.initial_energy_eV:.6f},{res.final_energy_eV:.6f},{res.max_energy_eV:.6f},{res.max_image_index},{res.n_images},{res.n_atoms}\n"
        )

    # Print concise human summary
    print("Au/Al(111) diffusion (fcc→fcc) via CI-NEB")
    print(f"Forward barrier E_d = {res.barrier_eV:.3f} eV (method: {res.barrier_method})")
    print(f"E_initial = {res.initial_energy_eV:.6f} eV, E_final = {res.final_energy_eV:.6f} eV, E_max = {res.max_energy_eV:.6f} eV")
    print(f"Images: {res.n_images}, atoms: {res.n_atoms}, peak image index: {res.max_image_index}")
    print(f"Energies per image (eV): {[round(x, 4) for x in energies]}")


if __name__ == "__main__":
    main()
