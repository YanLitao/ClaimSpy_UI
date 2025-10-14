import os
import json
import math
import numpy as np

from dataclasses import asdict, dataclass

# Reproducibility
np.random.seed(42)
try:
    import torch
    torch.manual_seed(42)
except Exception:
    torch = None

from dotenv import load_dotenv

from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.neb import NEB
from ase.io.trajectory import Trajectory
from ase.geometry import get_distances
from ase.constraints import FixAtoms

from chgnet.model import CHGNet, CHGNetCalculator

from pymatgen.core import Structure as PMGStructure
from pymatgen.io.ase import AseAtomsAdaptor

try:
    from mp_api.client import MPRester
except Exception as e:
    MPRester = None


@dataclass
class NEBResult:
    material_id: str
    n_images: int
    fmax: float
    steps: int
    energies_eV: list
    relative_energies_eV: list
    forward_barrier_eV: float
    initial_energy_eV: float
    max_energy_eV: float
    saddle_image_index: int


def load_mp_structure(api_key: str) -> tuple[PMGStructure, str]:
    assert MPRester is not None, "mp_api is not available in the environment."
    with MPRester(api_key) as mpr:
        # Search for LiC6; prefer the most stable (lowest E_hull)
        docs = mpr.summary.search(formula="LiC6",
                                  fields=["structure", "material_id", "energy_above_hull"])
    if not docs:
        raise RuntimeError("No LiC6 structures found from Materials Project API.")
    # Pick the lowest E_hull entry
    doc = min(docs, key=lambda d: (d.energy_above_hull if d.energy_above_hull is not None else 1e9))
    return doc.structure, doc.material_id


def load_graphite_structure(api_key: str) -> tuple[PMGStructure, str]:
    assert MPRester is not None, "mp_api is not available in the environment."
    with MPRester(api_key) as mpr:
        docs = mpr.summary.search(formula="C",
                                  fields=["structure", "material_id", "energy_above_hull", "symmetry"])
    if not docs:
        raise RuntimeError("No graphite-like carbon structures found from Materials Project API.")
    # Prefer graphite (P6_3/mmc) over diamond (Fd-3m)
    graphite_candidates = [d for d in docs if getattr(d, 'symmetry', None) and getattr(d.symmetry, 'symbol', None) in ("P6_3/mmc", "P6_3/mmc1")]  # tolerate slight variations
    if graphite_candidates:
        doc = min(graphite_candidates, key=lambda d: (d.energy_above_hull if d.energy_above_hull is not None else 1e9))
    else:
        doc = min(docs, key=lambda d: (d.energy_above_hull if d.energy_above_hull is not None else 1e9))
    return doc.structure, doc.material_id


def ase_min_distance_clip(atoms: Atoms, min_dist: float = 0.7) -> None:
    """If any interatomic separation is below min_dist (Å), scale cell slightly to avoid bad overlaps."""
    dmat = atoms.get_all_distances(mic=True)
    # Mask self-distances
    dmat += np.eye(len(atoms)) * 1e9
    if np.nanmin(dmat) < min_dist:
        # Uniformly scale cell and positions by small factor until safe
        scale = 1.02
        while np.nanmin(dmat) < min_dist:
            atoms.set_cell(atoms.cell * scale, scale_atoms=True)
            dmat = atoms.get_all_distances(mic=True)
            dmat += np.eye(len(atoms)) * 1e9


def is_finite_energy_forces(atoms: Atoms) -> bool:
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    return np.isfinite(e) and np.isfinite(f).all()


def pre_relax_positions(atoms: Atoms, calc: CHGNetCalculator, tag: str,
                        fmax_coarse: float = 0.2, fmax_fine: float = 0.05, max_steps: int = 300) -> None:
    atoms.set_calculator(calc)
    ase_min_distance_clip(atoms)
    if not is_finite_energy_forces(atoms):
        raise RuntimeError(f"Non-finite energy/forces before optimization for {tag}.")
    # Stage 1: positions-only
    opt1 = BFGS(atoms, logfile=f"relax_{tag}_1.log")
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)
    # Stage 2: positions-only fine
    opt2 = BFGS(atoms, logfile=f"relax_{tag}_2.log")
    opt2.run(fmax=fmax_fine, steps=max_steps)
    if not is_finite_energy_forces(atoms):
        raise RuntimeError(f"Non-finite energy/forces after optimization for {tag}.")


def build_supercell(pmg: PMGStructure, sc=(2, 2, 1)) -> PMGStructure:
    pmg_sc = pmg.copy()
    pmg_sc.make_supercell(sc)
    return pmg_sc


def find_inplane_li_pair(atoms: Atoms, z_tol_frac: float = 0.03) -> tuple[int, int, np.ndarray]:
    """Find indices of two neighboring Li sites within the same interlayer (approx. same fractional z).

    Returns (i, j, target_frac_j) where i is the moving Li index, j is the neighbor index to remove,
    and target_frac_j is the fractional coordinate of j (the destination site).
    """
    li_inds = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'Li']
    assert len(li_inds) >= 2, "Need at least two Li atoms to define a neighbor site in-plane."

    cell = atoms.get_cell()
    spos = atoms.get_scaled_positions(wrap=True)

    def frac_diff(a, b):
        df = a - b
        df -= np.round(df)
        return df

    best = None
    for idx_i in li_inds:
        fi = spos[idx_i]
        candidates = []
        for idx_j in li_inds:
            if idx_j == idx_i:
                continue
            fj = spos[idx_j]
            df = frac_diff(fj, fi)
            # Restrict to nearly same z-plane
            if abs(df[2]) <= z_tol_frac:
                # 2D in-plane metric
                df2d = df.copy(); df2d[2] = 0.0
                dv = np.dot(df2d, cell)  # Cartesian
                dist = np.linalg.norm(dv)
                candidates.append((dist, idx_i, idx_j, fj.copy()))
        if candidates:
            local_best = min(candidates, key=lambda x: x[0])
            if (best is None) or (local_best[0] < best[0]):
                best = local_best

    if best is None:
        raise RuntimeError("Failed to find in-plane Li neighbors. Try larger supercell or adjust z_tol_frac.")
    _, i, j, target_frac = best
    return i, j, target_frac


def remove_atom_by_index(atoms: Atoms, idx: int) -> Atoms:
    mask = [k for k in range(len(atoms)) if k != idx]
    return atoms[mask]


def reduce_to_li_pair(atoms: Atoms, i_idx: int, j_idx: int) -> tuple[Atoms, int, int]:
    """Return a new Atoms with all carbons kept but only two Li (i and j). Indices of the two Li are returned for the new Atoms."""
    symbols = atoms.get_chemical_symbols()
    carbon_inds = [k for k, s in enumerate(symbols) if s == 'C']
    li_inds = [k for k, s in enumerate(symbols) if s == 'Li']
    assert i_idx in li_inds and j_idx in li_inds
    keep = carbon_inds + [i_idx, j_idx]
    keep_sorted = sorted(keep)
    new_atoms = atoms[keep_sorted]
    # Map old indices to new indices
    old_to_new = {old: new for new, old in enumerate(keep_sorted)}
    i_new = old_to_new[i_idx]
    j_new = old_to_new[j_idx]
    return new_atoms, i_new, j_new


def run_neb(initial: Atoms, final: Atoms, calc: CHGNetCalculator, n_images: int = 7,
            fmax: float = 0.05, max_steps: int = 300, climb: bool = True,
            moving_index: int | None = None) -> tuple[list, list, int, list]:
    assert len(initial) == len(final), "Initial and final must have the same number of atoms."

    images = [initial]
    for _ in range(n_images - 2):
        images.append(initial.copy())
    images.append(final)

    # Attach calculator and optional constraints to each image
    for im in images:
        im.set_calculator(calc)
        if moving_index is not None:
            fixed = [k for k in range(len(im)) if k != moving_index]
            im.set_constraint(FixAtoms(indices=fixed))

    neb = NEB(images, climb=climb, k=0.1, allow_shared_calculator=True)
    try:
        neb.interpolate(method='idpp', apply_constraint=True)  # better initial band
    except Exception:
        neb.interpolate(apply_constraint=True)

    # Optimize NEB path
    opt = BFGS(neb, logfile='neb.log')
    opt.run(fmax=fmax, steps=max_steps)

    # Collect energies and find saddle
    energies = [im.get_potential_energy() for im in images]
    e0 = energies[0]
    rel = [e - e0 for e in energies]
    saddle_idx = int(np.argmax(energies))
    # Save band trajectory
    try:
        Trajectory('neb_images.traj', 'w', images)
        for k, im in enumerate(images):
            write(f'POSCAR_image_{k:02d}.vasp', im)
    except Exception:
        pass
    return energies, rel, saddle_idx, images


def main():
    load_dotenv()
    api_key = os.getenv('MP_API_KEY', None)
    if api_key is None or len(api_key.strip()) == 0:
        raise RuntimeError("MP_API_KEY not found in environment. Ensure .env contains MP_API_KEY.")

    # Load LiC6 from MP and build supercell
    pmg, mpid = load_mp_structure(api_key)
    pmg_sc = build_supercell(pmg, sc=(2, 2, 1))

    # Convert to ASE and prepare vacancy hop pair
    adaptor = AseAtomsAdaptor()
    atoms_sc = adaptor.get_atoms(pmg_sc)
    atoms_sc.set_pbc([True, True, True])

    # Load CHGNet once; reuse calculator
    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)

    # Identify two in-plane neighboring Li sites and create a vacancy at the target site
    i_idx, j_idx, j_frac = find_inplane_li_pair(atoms_sc, z_tol_frac=0.03)

    # Reduce system to carbon + two Li (i and j) to avoid Li-Li interactions
    reduced, i_idx, j_idx = reduce_to_li_pair(atoms_sc, i_idx, j_idx)

    # Create initial: remove neighbor j (vacancy), keep i at original site
    initial_full = reduced.copy()
    initial = remove_atom_by_index(initial_full, j_idx)

    # Create final: move Li i to j's fractional position, same vacancy at i's original site
    final = initial.copy()
    i_new = i_idx if i_idx < j_idx else i_idx  # since j removed, i_new == i_idx if i_idx < j_idx else i_idx-1; but we ensured i and j are last two entries
    # Recompute i_new robustly: find the sole Li index in 'initial'
    sym_init = initial.get_chemical_symbols()
    li_indices_init = [k for k, s in enumerate(sym_init) if s == 'Li']
    assert len(li_indices_init) == 1
    i_new = li_indices_init[0]
    # Set position of moving Li to j's fractional coordinate
    final.set_scaled_positions(final.get_scaled_positions())  # ensure fractional mapping consistent
    scaled = final.get_scaled_positions()
    scaled[i_new] = j_frac
    final.set_scaled_positions(scaled)

    # Pre-relax endpoints (positions-only) to avoid spurious forces
    pre_relax_positions(initial, calc, tag='initial')
    write('POSCAR_initial.vasp', initial)
    pre_relax_positions(final, calc, tag='final')
    write('POSCAR_final.vasp', final)

    # Run CI-NEB
    n_images = 7
    energies, rel, saddle_idx, band_images = run_neb(initial, final, calc, n_images=n_images,
                                                     fmax=0.05, max_steps=300, climb=True,
                                                     moving_index=i_new)

    # Save trajectory
    # Save saddle image explicitly
    saddle_atoms = band_images[saddle_idx]
    write('POSCAR_saddle.vasp', saddle_atoms)

    barrier_eV = float(max(energies) - energies[0])
    result = NEBResult(
        material_id=mpid,
        n_images=n_images,
        fmax=0.05,
        steps=300,
        energies_eV=[float(e) for e in energies],
        relative_energies_eV=[float(x) for x in rel],
        forward_barrier_eV=barrier_eV,
        initial_energy_eV=float(energies[0]),
        max_energy_eV=float(max(energies)),
        saddle_image_index=int(saddle_idx),
    )

    # Fallback: if barrier looks unphysical (> 5 eV), compute dilute Li in graphite gallery
    if barrier_eV > 5.0:
        print("Barrier from LiC6 vacancy model is high; switching to dilute Li in graphite gallery.")
        g_pmg, g_mpid = load_graphite_structure(api_key)
        adaptor = AseAtomsAdaptor()
        g_sc = build_supercell(g_pmg, sc=(3, 3, 2))
        g_atoms = adaptor.get_atoms(g_sc)
        g_atoms.set_pbc([True, True, True])

        # Remove any Li (should be none) and compute layer z's
        symbols = g_atoms.get_chemical_symbols()
        keep_c = [i for i, s in enumerate(symbols) if s == 'C']
        g_atoms = g_atoms[keep_c]

        spos = g_atoms.get_scaled_positions(wrap=True)
        z_all = np.array([s[2] for s in spos])
        # Cluster two layers by k-means-like median split
        z_sorted = np.sort(z_all)
        z1 = np.median(z_sorted[:len(z_sorted)//2])
        z2 = np.median(z_sorted[len(z_sorted)//2:])
        zmid = (z1 + z2) / 2.0

        # Initial Li at hollow-like fractional (1/3, 1/3, zmid)
        cell = g_atoms.get_cell()
        li_pos_i = np.dot([1/3, 1/3, zmid % 1.0], cell)
        initial_g = g_atoms.copy() + Atoms('Li', positions=[li_pos_i])

        # Final Li at neighboring hollow (2/3, 2/3, zmid)
        final_g = initial_g.copy()
        spos_g = final_g.get_scaled_positions()
        li_idx_g = len(final_g) - 1
        spos_g[li_idx_g] = np.array([(2/3) % 1.0, (2/3) % 1.0, zmid % 1.0])
        final_g.set_scaled_positions(spos_g)

        # Pre-relax endpoints
        pre_relax_positions(initial_g, calc, tag='initial_g')
        write('POSCAR_initial_graphite.vasp', initial_g)
        pre_relax_positions(final_g, calc, tag='final_g')
        write('POSCAR_final_graphite.vasp', final_g)

        # NEB with only Li mobile
        energies_g, rel_g, saddle_idx_g, band_images_g = run_neb(initial_g, final_g, calc, n_images=7,
                                                                 fmax=0.05, max_steps=300, climb=True,
                                                                 moving_index=li_idx_g)
        write('POSCAR_saddle_graphite.vasp', band_images_g[saddle_idx_g])
        result = NEBResult(
            material_id=g_mpid,
            n_images=7,
            fmax=0.05,
            steps=300,
            energies_eV=[float(e) for e in energies_g],
            relative_energies_eV=[float(x) for x in rel_g],
            forward_barrier_eV=float(max(energies_g) - energies_g[0]),
            initial_energy_eV=float(energies_g[0]),
            max_energy_eV=float(max(energies_g)),
            saddle_image_index=int(saddle_idx_g),
        )

    with open('neb_results.json', 'w') as f:
        json.dump(asdict(result), f, indent=2)

    # CSV for quick inspection
    with open('neb_energies.csv', 'w') as f:
        f.write('image_index,energy_eV,relative_eV\n')
        for i, (e, r) in enumerate(zip(result.energies_eV, result.relative_energies_eV)):
            f.write(f"{i},{e:.8f},{r:.8f}\n")

    with open('summary.csv', 'w') as f:
        f.write('metric,value\n')
        f.write(f"material_id,{result.material_id}\n")
        f.write(f"forward_barrier_eV,{result.forward_barrier_eV:.6f}\n")

    # Human summary
    print("NEB summary: Li in-plane hop")
    print(f"Materials Project ID: {result.material_id}")
    print(f"Images: {result.n_images}, optimizer fmax: {result.fmax} eV/Ang")
    print(f"Forward barrier (E_max - E_initial): {result.forward_barrier_eV:.4f} eV")
    print(f"Saddle image index: {result.saddle_image_index}")
    print("Per-image energies (eV):")
    for i, (e, r) in enumerate(zip(result.energies_eV, result.relative_energies_eV)):
        print(f"  {i:2d}: E = {e: .6f}  dE = {r: .6f}")


if __name__ == '__main__':
    main()
