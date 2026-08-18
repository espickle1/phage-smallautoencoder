"""Top-5 SAE features per HP3 ROI residue, with descriptions from the biohub
API.

ROI set = the 8 rows already in ``analysis/structures/masked_marginals_hp3.csv``
(the established source of truth for HP3 ROIs): tcf's 3 substitutions in WT
context, gp8's WT flank + 3 inserted residues + evolved flank. Activation
lookup is fully local (``analysis/features/<protein>_<variant>.npz``, key
``features``); descriptions/categories are one GET per unique feature id to
``https://biohub.ai/esm/protein/api/v1alpha1/features/<id>`` -- public, no
token needed (confirmed 2026-07-21).

Usage:
    python src/hp3_feature_lookup.py
Writes analysis/features/hp3_roi_features_top5.csv.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from hp3_feature_cache import get_feature_info

ROOT = Path(__file__).resolve().parent.parent
FEATURES = ROOT / "analysis" / "features"
MASKED_MARGINALS = ROOT / "analysis" / "structures" / "masked_marginals_hp3.csv"
OUT = FEATURES / "hp3_roi_features_top5.csv"
TOP_K = 5

PROTEIN_NAME = {"gp8": "gp8", "tcf": "tail_collar_fiber"}


def top_features(protein_key: str, variant: str, position: int, top_k: int = TOP_K):
    npz = np.load(FEATURES / f"{protein_key}_{variant}.npz")
    acts = npz["features"][position - 1]  # 1-indexed position -> 0-indexed row
    idx = np.argsort(acts)[::-1]
    idx = idx[acts[idx] > 0][:top_k]
    return [(int(i), float(acts[i])) for i in idx]


def main() -> None:
    mm = pd.read_csv(MASKED_MARGINALS)
    rows = []
    for r in mm.itertuples():
        feats = top_features(r.protein, r.context_variant, int(r.position))
        for rank, (feature_id, activation) in enumerate(feats, 1):
            info = get_feature_info(feature_id)
            rows.append(dict(
                protein=PROTEIN_NAME[r.protein], variant=r.context_variant,
                position=int(r.position), residue=r.original, kind=r.kind,
                rank=rank, feature_id=feature_id, activation=activation,
                label=info.get("label", ""), category=info.get("category", "Uncategorized"),
            ))
        print(f"{r.protein}/{r.context_variant} {r.original}{r.position}: "
              f"{len(feats)} features fetched")

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(ROOT)} ({len(df)} rows, "
          f"{df['feature_id'].nunique()} unique features)")


if __name__ == "__main__":
    main()
