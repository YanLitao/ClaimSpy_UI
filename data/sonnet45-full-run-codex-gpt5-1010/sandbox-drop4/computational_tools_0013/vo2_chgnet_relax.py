import os
import json
import csv
import math
import random
from dataclasses import asdict, dataclass
from typing import List, Dict, Any

import numpy as np

# Load API keys
from dotenv import load_dotenv

# Structure handling
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.io.ase import AseAtomsAdaptor

# CHGNet / ASE
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.atoms import Atoms


def set_seeds(seed: int = 42):
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass


def load_mp_structure_vo2() -> Structure:
    """Fetch a low-temperature VO2 structure (monoclinic if available) from Materials Project.

    Prefers stable entries and lowest energy above hull. Falls back to any VO2 if filtering fails.
    Requires MP_API_KEY in environment.
    """
    api_key = os.getenv("MP_API_KEY")
    if not api_key:
        raise RuntimeError("MP_API_KEY not found. Ensure .env contains MP_API_KEY and is loaded.")

    # Try new mp-api client first
    try:
        from mp_api.client import MPRester as MPClient

        with MPClient(api_key=api_key) as mpr:
            docs = mpr.materials.summary.search(
                formula="VO2",
                is_stable=True,
                fields=[
                    "material_id",
                    "structure",
                    "energy_above_hull",
                    "symmetry",
                ],
            )

        if not docs:
            # Fallback: drop is_stable filter
            with MPClient(api_key=api_key) as mpr:
                docs = mpr.materials.summary.search(
                    formula="VO2",
                    fields=["material_id", "structure", "energy_above_hull", "symmetry"],
                )

        if not docs:
            raise RuntimeError("No VO2 documents returned from MP API.")

        # Prefer P21/c monoclinic (M1) if present; else lowest energy_above_hull
        def is_m1(doc):
            try:
                sym = getattr(doc, "symmetry", None)
                if sym is None and isinstance(doc, dict):
                    sym = doc.get("symmetry", {})
                symbol = getattr(sym, "symbol", None) if sym is not None else None
                if symbol is None and isinstance(sym, dict):
                    symbol = sym.get("symbol")
                return symbol in ("P2_1/c", "P21/c", "P 1 2/c 1", "P 1 2/c 1 (unique b)")
            except Exception:
                return False

        m1_docs = [d for d in docs if is_m1(d)]
        if m1_docs:
            # choose one with minimum energy_above_hull
            chosen = min(
                m1_docs, key=lambda d: getattr(d, "energy_above_hull", getattr(d, "energy_above_hull", 0.0))
            )
        else:
            chosen = min(
                docs, key=lambda d: getattr(d, "energy_above_hull", getattr(d, "energy_above_hull", 0.0))
            )

        structure = chosen.structure if hasattr(chosen, "structure") else chosen["structure"]
        if not isinstance(structure, Structure):
            # mp-api usually returns Pymatgen Structure; if not, attempt construction
            structure = Structure.from_dict(structure)
        return structure

    except Exception:
        # Fallback to legacy pymatgen MPRester
        from pymatgen.ext.matproj import MPRester as LegacyMR

        with LegacyMR(api_key) as mr:
            # Try to fetch stable VO2 entries; request necessary fields
            results = mr.query(
                {"formula_pretty": "VO2", "deprecated": False},
                ["material_id", "structure", "e_above_hull", "spacegroup.symbol", "is_stable"],
            )
        if not results:
            raise RuntimeError("No VO2 results from legacy MPRester.")

        # Filter stable first
        stable = [r for r in results if r.get("is_stable")]
        pool = stable if stable else results

        def is_m1_res(r):
            sym = r.get("spacegroup.symbol") or r.get("spacegroup", {}).get("symbol")
            return sym in ("P2_1/c", "P21/c", "P 1 2/c 1")

        m1_pool = [r for r in pool if is_m1_res(r)]
        if m1_pool:
            chosen = min(m1_pool, key=lambda r: r.get("e_above_hull", 0.0))
        else:
            chosen = min(pool, key=lambda r: r.get("e_above_hull", 0.0))

        structure = chosen["structure"]
        if not isinstance(structure, Structure):
            structure = Structure.from_dict(structure)
        return structure


def ensure_min_distances(atoms: Atoms, min_dist: float = 0.7) -> None:
    """If any interatomic distance is < min_dist Å, scale the cell slightly to avoid clashes."""
    dmat = atoms.get_all_distances(mic=True)
    np.fill_diagonal(dmat, np.inf)
    dmin = float(np.min(dmat))
    if dmin < min_dist:
        scale = min_dist / max(dmin, 1e-6)
        # Add margin
        scale = max(scale, 1.05)
        cell = atoms.cell.array
        atoms.set_cell(cell * scale, scale_atoms=True)


def check_finite(atoms: Atoms, stage: str) -> None:
    e = atoms.get_potential_energy()
    f = atoms.get_forces()
    if not np.isfinite(e) or not np.isfinite(f).all():
        raise RuntimeError(f"Non-finite energy/forces {stage}; check initial structure.")


def relax_vo2(structure: Structure, workdir: str) -> Atoms:
    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(structure)

    # Load CHGNet once
    from chgnet.model import CHGNet, CHGNetCalculator

    chgnet = CHGNet.load()
    calc = CHGNetCalculator(chgnet)
    atoms.set_calculator(calc)

    ensure_min_distances(atoms, min_dist=0.7)

    # Pre-check
    check_finite(atoms, stage="before relaxation")

    # Stage 1: positions only
    fmax_coarse = 0.2
    max_steps = 300
    log1 = os.path.join(workdir, "relax1.log")
    opt1 = BFGS(atoms, logfile=log1)
    opt1.run(fmax=fmax_coarse, steps=max_steps // 2)

    check_finite(atoms, stage="after stage 1")

    # Stage 2: hydrostatic cell relax
    fmax_fine = 0.05
    log2 = os.path.join(workdir, "relax2.log")
    ucf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt2 = BFGS(ucf, logfile=log2)
    opt2.run(fmax=fmax_fine, steps=max_steps)

    check_finite(atoms, stage="after stage 2")

    return atoms


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    """Return angle between vectors in degrees."""
    a = np.dot(v1, v2)
    b = np.linalg.norm(v1) * np.linalg.norm(v2)
    if b == 0:
        return float("nan")
    c = np.clip(a / b, -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))


@dataclass
class OctahedronReport:
    v_index: int
    v_frac_coords: List[float]
    neighbor_o_indices: List[int]
    bond_lengths: List[float]
    bond_length_min: float
    bond_length_max: float
    bond_length_mean: float
    bond_length_std: float
    bond_length_range: float
    n_angles: int
    angle_min: float
    angle_max: float
    angle_mean: float
    angle_std: float
    angle_max_dev_from_90_180: float
    distorted: bool


def analyze_vo6(structure: Structure) -> Dict[str, Any]:
    species = [str(sp) for sp in structure.species]
    v_indices = [i for i, s in enumerate(species) if s.startswith("V")]
    o_indices = [i for i, s in enumerate(species) if s.startswith("O")]
    if not v_indices or len(o_indices) < 6:
        raise RuntimeError("Structure does not contain expected V and O counts for VO2.")

    # Compute V-O distance matrix under PBC using fractional coords
    dmat = structure.lattice.get_all_distances(
        structure.frac_coords[v_indices], structure.frac_coords[o_indices]
    )

    reports: List[OctahedronReport] = []

    length_range_threshold = 0.05  # Å
    angle_dev_threshold = 2.0  # degrees

    for iv_local, v_idx in enumerate(v_indices):
        row = dmat[iv_local]
        local_sorted = np.argsort(row)[:6]
        chosen_o_local = list(map(int, local_sorted))
        neighbor_o_indices = [o_indices[j] for j in chosen_o_local]

        # Use PMG to get consistent vectors with periodic images
        vecs = []
        bond_lengths = []
        for o_idx in neighbor_o_indices:
            dist, image = structure[v_idx].distance_and_image(structure[o_idx])
            dv_frac = structure.frac_coords[o_idx] + np.array(image) - structure.frac_coords[v_idx]
            dv_cart = structure.lattice.get_cartesian_coords(dv_frac)
            vecs.append(dv_cart)
            bond_lengths.append(float(dist))

        # Angles: all unique pairs
        angles = []
        for a in range(6):
            for b in range(a + 1, 6):
                ang = angle_between(vecs[a], vecs[b])
                angles.append(ang)

        # Stats
        bl = np.array(bond_lengths)
        ang = np.array(angles)

        # Deviation relative to nearest of 90 or 180
        dev_90 = np.abs(ang - 90.0)
        dev_180 = np.abs(ang - 180.0)
        dev = np.minimum(dev_90, dev_180)

        bond_length_min = float(np.min(bl))
        bond_length_max = float(np.max(bl))
        bond_length_mean = float(np.mean(bl))
        bond_length_std = float(np.std(bl))
        bond_length_range = float(bond_length_max - bond_length_min)

        angle_min = float(np.min(ang))
        angle_max = float(np.max(ang))
        angle_mean = float(np.mean(ang))
        angle_std = float(np.std(ang))
        angle_max_dev = float(np.max(dev))

        distorted = (bond_length_range >= length_range_threshold) or (angle_max_dev >= angle_dev_threshold)

        reports.append(
            OctahedronReport(
                v_index=int(v_idx),
                v_frac_coords=list(map(float, structure.frac_coords[v_idx])),
                neighbor_o_indices=list(map(int, neighbor_o_indices)),
                bond_lengths=list(map(float, bond_lengths)),
                bond_length_min=bond_length_min,
                bond_length_max=bond_length_max,
                bond_length_mean=bond_length_mean,
                bond_length_std=bond_length_std,
                bond_length_range=bond_length_range,
                n_angles=len(angles),
                angle_min=angle_min,
                angle_max=angle_max,
                angle_mean=angle_mean,
                angle_std=angle_std,
                angle_max_dev_from_90_180=angle_max_dev,
                distorted=distorted,
            )
        )

    overall_distorted = any(r.distorted for r in reports)
    return {
        "n_v": len(v_indices),
        "length_range_threshold_A": length_range_threshold,
        "angle_dev_threshold_deg": angle_dev_threshold,
        "overall_distorted": overall_distorted,
        "reports": [asdict(r) for r in reports],
    }


def write_outputs(
    workdir: str,
    initial: Structure,
    final_atoms: Atoms,
    analysis: Dict[str, Any],
):
    adaptor = AseAtomsAdaptor()
    final_struct = adaptor.get_structure(final_atoms)

    # Files
    poscar_init = os.path.join(workdir, "POSCAR.vo2_initial")
    poscar_final = os.path.join(workdir, "POSCAR.vo2_final")
    cif_final = os.path.join(workdir, "vo2_final.cif")
    json_out = os.path.join(workdir, "vo2_relax_results.json")
    csv_out = os.path.join(workdir, "vo2_relax_summary.csv")

    initial.to(fmt="poscar", filename=poscar_init)
    final_struct.to(fmt="poscar", filename=poscar_final)
    final_struct.to(fmt="cif", filename=cif_final)

    # JSON
    with open(json_out, "w") as f:
        json.dump(analysis, f, indent=2)

    # CSV summary per V
    headers = [
        "v_index",
        "bond_length_min_A",
        "bond_length_max_A",
        "bond_length_mean_A",
        "bond_length_std_A",
        "bond_length_range_A",
        "angle_min_deg",
        "angle_max_deg",
        "angle_mean_deg",
        "angle_std_deg",
        "angle_max_dev_from_90_180_deg",
        "distorted",
    ]
    with open(csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in analysis["reports"]:
            w.writerow(
                [
                    r["v_index"],
                    f"{r['bond_length_min']:.4f}",
                    f"{r['bond_length_max']:.4f}",
                    f"{r['bond_length_mean']:.4f}",
                    f"{r['bond_length_std']:.4f}",
                    f"{r['bond_length_range']:.4f}",
                    f"{r['angle_min']:.2f}",
                    f"{r['angle_max']:.2f}",
                    f"{r['angle_mean']:.2f}",
                    f"{r['angle_std']:.2f}",
                    f"{r['angle_max_dev_from_90_180']:.2f}",
                    r["distorted"],
                ]
            )

    return {
        "poscar_initial": poscar_init,
        "poscar_final": poscar_final,
        "cif_final": cif_final,
        "json": json_out,
        "csv": csv_out,
        "relax_logs": [
            os.path.join(workdir, "relax1.log"),
            os.path.join(workdir, "relax2.log"),
        ],
    }


def main():
    set_seeds(42)
    load_dotenv(dotenv_path=os.path.join(os.getcwd(), ".env"))
    workdir = os.getcwd()

    print("Fetching VO2 structure from MP...")
    initial_struct = load_mp_structure_vo2()
    formula = initial_struct.composition.reduced_formula
    print(f"Fetched structure: {formula}, n_sites={len(initial_struct)}")
    print(f"Initial cell (Å): a,b,c = {initial_struct.lattice.a:.3f}, {initial_struct.lattice.b:.3f}, {initial_struct.lattice.c:.3f}; alpha,beta,gamma = {initial_struct.lattice.alpha:.2f}, {initial_struct.lattice.beta:.2f}, {initial_struct.lattice.gamma:.2f}")

    print("Starting CHGNet relaxation (positions then hydrostatic cell)...")
    final_atoms = relax_vo2(initial_struct, workdir)
    final_energy = final_atoms.get_potential_energy()
    print(f"Final potential energy (eV): {final_energy:.6f}")

    print("Analyzing VO6 octahedra...")
    analysis = analyze_vo6(AseAtomsAdaptor().get_structure(final_atoms))

    # Human-readable summary
    print("--- VO6 Octahedra Summary (per V) ---")
    print("v_index, range_A, std_A, max_angle_dev_deg, distorted")
    for r in analysis["reports"]:
        print(
            f"{r['v_index']}, {r['bond_length_range']:.4f}, {r['bond_length_std']:.4f}, {r['angle_max_dev_from_90_180']:.2f}, {r['distorted']}"
        )
    print(f"Overall distorted: {analysis['overall_distorted']}")

    files = write_outputs(workdir, initial_struct, final_atoms, analysis)

    # Final concise conclusion
    if analysis["overall_distorted"]:
        conclusion = (
            "Conclusion: VO6 octahedra are distorted (non-uniform V–O bonds and/or O–V–O angles deviate from ideal)."
        )
    else:
        conclusion = (
            "Conclusion: VO6 octahedra are not detectably distorted within thresholds (<=0.05 Å range, <=2° dev)."
        )
    print(conclusion)

    # Emit a minimal machine-readable pointer to outputs
    manifest = {
        "poscar_initial": files["poscar_initial"],
        "poscar_final": files["poscar_final"],
        "cif_final": files["cif_final"],
        "json": files["json"],
        "csv": files["csv"],
        "relax_logs": files["relax_logs"],
    }
    with open(os.path.join(workdir, "vo2_file_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)


if __name__ == "__main__":
    main()
