"""Minimum heavy-atom distance from each tail_collar_fiber ROI residue to the
OTHER two chains of its homotrimer, in both WT and evolved AF3 trimer models.

WHY: tail_collar_fiber's real biological unit is a homotrimer (3 identical
chains), not the isolated monomer used elsewhere in this HP3 pipeline. The
mutations (Y444H, K464R, Y465H) only matter structurally if they sit near an
inter-chain interface.

Structures (verified 2026-07-21): both are true homotrimers, 3 identical
516-residue chains (A/B/C), auth_seq_id 1..516 contiguous, no offset. The
evolved trimer was predicted AFTER the Y444H FASTA correction (all 3 chains
carry H444/R464/H465), unlike the earlier monomer evolved structure -- so
this analysis (and the trimer figures built from it) has no stale-residue
caveat.

Usage:
    python src/tcf_interface_contacts.py
Writes analysis/structures/tcf_interface_contacts.csv.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
STRUCT_DIR = ROOT / "data" / "structures"
OUT = ROOT / "analysis" / "structures" / "tcf_interface_contacts.csv"

TRIMERS = {
    "wt": STRUCT_DIR / "fold_yp_010228801_tail_collar_fiber_3mer_wt_model_0.cif",
    "evolved": STRUCT_DIR / "fold_yp_010228801_tail_collar_fiber_3mer_hp3e_model_0.cif",
}
CHAINS = ["A", "B", "C"]
ROI_POSITIONS = [444, 464, 465]
ROI_RESIDUE = {"wt": {444: "Y", 464: "K", 465: "Y"},
               "evolved": {444: "H", 464: "R", 465: "H"}}


def load_atoms(path: Path) -> dict[tuple[str, int], np.ndarray]:
    """(chain, position) -> [n_atoms, 3] heavy-atom coordinates."""
    atoms: dict[tuple[str, int], list] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        chain, pos = f[16], int(f[15])
        xyz = (float(f[10]), float(f[11]), float(f[12]))
        atoms.setdefault((chain, pos), []).append(xyz)
    return {k: np.array(v) for k, v in atoms.items()}


def min_dist(a: np.ndarray, b: np.ndarray) -> float:
    diff = a[:, None, :] - b[None, :, :]
    return float(np.sqrt((diff ** 2).sum(-1)).min())


def main() -> None:
    rows = []
    for variant, path in TRIMERS.items():
        atoms = load_atoms(path)
        # all heavy-atom coords per chain, for nearest-atom-anywhere-on-that-chain
        # lookups (the interface partner isn't necessarily the same residue
        # number on the other chain).
        chain_atoms = {c: np.vstack([v for (c2, _p2), v in atoms.items() if c2 == c])
                       for c in CHAINS}
        for pos in ROI_POSITIONS:
            for chain in CHAINS:
                self_atoms = atoms[(chain, pos)]
                dists = {oc: min_dist(self_atoms, chain_atoms[oc])
                         for oc in CHAINS if oc != chain}
                nearest_chain = min(dists, key=dists.get)
                rows.append(dict(
                    variant=variant, position=pos, residue=ROI_RESIDUE[variant][pos],
                    chain=chain, nearest_chain=nearest_chain,
                    min_distance_A=round(dists[nearest_chain], 2),
                    **{f"dist_to_{oc}_A": round(d, 2) for oc, d in dists.items()},
                ))
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(df.to_string(index=False))
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
