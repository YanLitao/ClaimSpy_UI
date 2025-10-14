"""
Zr–Sb A15 (Zr3Sb) stability check via CALPHAD

Follows CALPHAD.md best practices, with a fast fail if the local
TDB inventory lacks the required elements.

Outputs a concise JSON summary for downstream consumption and prints
an explainer suitable for a research log.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List


def main() -> None:
    # --- Cached minimal info from CALPHAD.md (read once elsewhere) ---
    # Local TDBs reside under the repo root in a folder named 'TnETDBDB'.
    cwd = Path(__file__).resolve().parent
    # Find repo root by climbing until 'TnETDBDB' exists as a child
    repo_root = None
    cur = cwd
    for _ in range(6):  # limit climb depth
        candidate = cur / 'TnETDBDB'
        if candidate.is_dir():
            repo_root = cur
            break
        cur = cur.parent
    if repo_root is None:
        # Fall back to two levels up (expected in typical layout)
        repo_root = cwd.parents[1]
    tdb_dir = repo_root / 'TnETDBDB'

    claim = (
        "Is A15 Zr3Sb (25 at.% Sb) thermodynamically stable at 1 atm? "
        "Report convex hull distance (meV/atom), competing phases, and any T stability window."
    )

    out: Dict[str, Any] = {
        "claim": claim,
        "system": "Zr–Sb",
        "composition_target": {"ZR": 0.75, "SB": 0.25},
        "pressure_Pa": 101325,
        "provenance": {
            "tdb_dir": str(tdb_dir),
            "tdb_candidates": [],
            "tdb_selected": None,
        },
        "result": None,
        "error": None,
    }

    if not tdb_dir.is_dir():
        msg = (
            f"Local TDB directory not found at {tdb_dir}. Cannot perform CALPHAD "
            f"for Zr–Sb under the repository constraints."
        )
        out["error"] = {"type": "MissingTDBDir", "message": msg}
        _emit(out)
        return

    # Enumerate local TDBs (case-insensitive .tdb)
    tdb_paths = sorted([p for p in tdb_dir.iterdir() if p.suffix.lower() == ".tdb"])
    out["provenance"]["tdb_candidates"] = [p.name for p in tdb_paths]

    if not tdb_paths:
        msg = (
            f"No .tdb files found in {tdb_dir}. Unable to evaluate Zr–Sb stability." 
        )
        out["error"] = {"type": "NoTDBFiles", "message": msg}
        _emit(out)
        return

    # Try to load each TDB and check for required elements
    # We defer importing pycalphad until we know we have a candidate, to avoid slow imports needlessly.
    required = {"ZR", "SB", "VA"}
    selected = None
    try:
        from pycalphad import Database
    except Exception as e:  # pragma: no cover
        out["error"] = {"type": "ImportError", "message": f"pycalphad import failed: {e}"}
        _emit(out)
        return

    for tdb in tdb_paths:
        try:
            db = Database(str(tdb))
            present = {str(el).upper() for el in getattr(db, 'elements', [])}
            if required.issubset(present):
                selected = tdb
                break
        except Exception:
            # Skip unreadable/malformed candidates
            continue

    if selected is None:
        msg = (
            "Available local TDBs do not provide the required element set "
            "{Zr, Sb, VA}. The present database(s) cannot represent the Zr–Sb binary. "
            "Per CALPHAD.md, no external TDBs may be used, so this task cannot be computed here."
        )
        out["error"] = {"type": "ElementsNotCovered", "message": msg}
        _emit(out)
        return

    # If we ever reach here, we would proceed with full equilibrium/hull analysis.
    out["provenance"]["tdb_selected"] = selected.name
    out["result"] = {
        "note": (
            "Placeholder: a suitable TDB was found. In the current repository snapshot, "
            "this branch is not expected to execute because Sb is absent from COST507R."
        )
    }
    _emit(out)


def _emit(payload: Dict[str, Any]) -> None:
    """Write results JSON and print a concise human summary."""
    outpath = Path("zr_sb_a15_results.json")
    outpath.write_text(json.dumps(payload, indent=2))

    if payload.get("error"):
        err = payload["error"]["message"]
        # Human-facing summary
        print("CALPHAD run aborted — local database limitation.")
        print(err)
        print("Candidates:", ", ".join(payload["provenance"].get("tdb_candidates", [])))
        print(f"Wrote: {outpath}")
        # Exit with non-zero for clarity in scripted contexts
        sys.exit(2)
    else:
        print("CALPHAD preliminary selection complete.")
        print(json.dumps(payload, indent=2))
        print(f"Wrote: {outpath}")


if __name__ == "__main__":
    main()
