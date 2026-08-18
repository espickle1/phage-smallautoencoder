"""Whole-sequence top-3 SAE feature scan for HP3 gp8 / tail_collar_fiber --
run one protein at a time, per the 'fully separate analyses' decision.

For every residue in both the WT and evolved sequence, records the top-3
active SAE features (by activation) directly from the local
``analysis/features/<protein>_<variant>.npz`` arrays (no network). The set of
DISTINCT feature ids encountered is then looked up once each via
``hp3_feature_cache`` (disk-persisted -- ids already seen from the ROI/cluster
work are free; only genuinely new ids hit the network).

Outputs
-------
    analysis/features/hp3_whole_protein_features_top3_<protein>.csv
        long format: protein, variant, position, residue, rank, feature_id,
        activation, label, category -- one row per (residue, rank).
    analysis/features/hp3_feature_catalog.csv
        deduplicated feature_id -> label/category/summary, shared reference
        across the whole project (grows as more HP3 scripts run).
    analysis/figures/feature_track_hp3_<protein>.png
        rank-1 (dominant) feature's category along the whole sequence, one
        row per variant, ROI positions marked -- the feature-level analog of
        entropy_track.py.

Also prints every position where WT and evolved DISAGREE on their rank-1
feature (using data/gp8_alignment.csv for gp8's post-insertion position map;
tail_collar_fiber has no indel, so it's a direct 1:1 comparison) -- surfaces
mutation-adjacent effects beyond the pre-annotated ROI.

Usage:
    python src/hp3_whole_protein_features.py --protein gp8
    python src/hp3_whole_protein_features.py --protein tcf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from hp3_feature_cache import get_feature_info

ROOT = Path(__file__).resolve().parent.parent
FEATURES = ROOT / "analysis" / "features"
FIGURES = ROOT / "analysis" / "figures"
STRUCTURES = ROOT / "analysis" / "structures"
SEQ_DIR = ROOT / "data" / "sequences"
GP8_ALIGNMENT = ROOT / "data" / "gp8_alignment.csv"

TOP_K = 3
VARIANTS = ["wt", "evolved"]
SEQ_FILES = {"gp8": SEQ_DIR / "baseplate_gp8_pair.fasta",
             "tcf": SEQ_DIR / "tail_collar_fiber.fasta"}
PROTEIN_NAME = {"gp8": "gp8", "tcf": "tail_collar_fiber"}

CATEGORY_ORDER = [
    "Domain", "Disorder", "Interaction site", "Catalytic function",
    "Compositional bias", "Ligand-binding site", "Structural motif",
    "Post-translational modification", "Sequence motif", "Other", "Uncategorized",
]
CATEGORY_COLOR = {
    "Domain": "#2a78d6", "Disorder": "#1baf7a", "Interaction site": "#eda100",
    "Catalytic function": "#008300", "Compositional bias": "#4a3aa7",
    "Ligand-binding site": "#e34948", "Structural motif": "#e87ba4",
    "Post-translational modification": "#00a3a3", "Sequence motif": "#8b5a2b",
    "Other": "#a35200", "Uncategorized": "#898781",
}
INK, MUTED = "#0b0b0b", "#898781"


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header, chunks = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def _roi_keys(protein_key: str) -> set[tuple[str, int]]:
    mm = pd.read_csv(STRUCTURES / "masked_marginals_hp3.csv")
    mm = mm[mm["protein"] == protein_key]
    if protein_key == "tcf":
        positions = set(mm["position"])
        return {(v, p) for v in VARIANTS for p in positions}
    return set(zip(mm["context_variant"], mm["position"]))


def top_k_features(row: np.ndarray, k: int = TOP_K) -> list[tuple[int, float]]:
    idx = np.argsort(row)[::-1]
    idx = idx[row[idx] > 0][:k]
    return [(int(i), float(row[i])) for i in idx]


def scan(protein_key: str) -> pd.DataFrame:
    records = _parse_fasta(SEQ_FILES[protein_key])
    seqs = {"wt": records[0][1], "evolved": records[1][1]}

    rows = []
    unique_ids: set[int] = set()
    for variant in VARIANTS:
        X = np.load(FEATURES / f"{protein_key}_{variant}.npz")["features"]
        seq = seqs[variant]
        for pos in range(1, len(seq) + 1):
            for rank, (fid, act) in enumerate(top_k_features(X[pos - 1]), 1):
                unique_ids.add(fid)
                rows.append(dict(protein=protein_key, variant=variant, position=pos,
                                  residue=seq[pos - 1], rank=rank, feature_id=fid, activation=act))
    print(f"[{protein_key}] {len(rows)} (residue, rank) rows; {len(unique_ids)} unique features")

    info_map = {}
    for i, fid in enumerate(sorted(unique_ids), 1):
        info_map[fid] = get_feature_info(fid)
        if i % 200 == 0 or i == len(unique_ids):
            print(f"  fetched descriptions: {i}/{len(unique_ids)}")

    df = pd.DataFrame(rows)
    df["label"] = df["feature_id"].map(lambda f: info_map[f].get("label", ""))
    df["category"] = df["feature_id"].map(
        lambda f: info_map[f].get("category") if info_map[f].get("category") in CATEGORY_COLOR
        else "Uncategorized")
    out = FEATURES / f"hp3_whole_protein_features_top3_{protein_key}.csv"
    df.to_csv(out, index=False)
    print(f"  wrote {out.relative_to(ROOT)}")
    return df


def update_catalog(df: pd.DataFrame) -> None:
    cat_path = FEATURES / "hp3_feature_catalog.csv"
    new = (df[["feature_id", "label", "category"]].drop_duplicates("feature_id")
           .sort_values("feature_id"))
    if cat_path.exists():
        existing = pd.read_csv(cat_path)
        combined = (pd.concat([existing, new]).drop_duplicates("feature_id")
                    .sort_values("feature_id"))
    else:
        combined = new
    combined.to_csv(cat_path, index=False)
    print(f"  wrote {cat_path.relative_to(ROOT)} ({len(combined)} unique features total)")


def flag_variant_disagreements(protein_key: str, df: pd.DataFrame) -> None:
    rank1 = df[df["rank"] == 1].set_index(["variant", "position"])["feature_id"]
    roi_keys = _roi_keys(protein_key)

    if protein_key == "tcf":
        pairs = [(p, p) for p in sorted(df.loc[df.variant == "wt", "position"].unique())]
    else:
        align = pd.read_csv(GP8_ALIGNMENT).dropna(subset=["wt_position"])
        pairs = [(int(r.wt_position), int(r.evolved_position)) for r in align.itertuples()]

    print(f"\n  WT vs evolved rank-1 feature disagreements ({protein_key}):")
    n_diff = 0
    for wt_pos, ev_pos in pairs:
        f_wt = rank1.get(("wt", wt_pos))
        f_ev = rank1.get(("evolved", ev_pos))
        if f_wt is not None and f_ev is not None and f_wt != f_ev:
            n_diff += 1
            tag = "  [ROI]" if ("wt", wt_pos) in roi_keys or ("evolved", ev_pos) in roi_keys else ""
            print(f"    wt pos {wt_pos} (feature {f_wt}) vs evolved pos {ev_pos} "
                  f"(feature {f_ev}){tag}")
    print(f"  {n_diff}/{len(pairs)} aligned positions disagree on rank-1 feature")


def plot_track(protein_key: str, df: pd.DataFrame) -> None:
    roi_keys = _roi_keys(protein_key)
    rank1 = df[df["rank"] == 1]

    fig, axes = plt.subplots(2, 1, figsize=(14, 3.2), sharex=False)
    for ax, variant in zip(axes, VARIANTS):
        sub = rank1[rank1.variant == variant].sort_values("position")
        colors = [CATEGORY_COLOR.get(c, CATEGORY_COLOR["Uncategorized"]) for c in sub["category"]]
        ax.bar(sub["position"], 1, width=1.0, color=colors, linewidth=0)
        roi_pos = [p for v, p in roi_keys if v == variant]
        for p in roi_pos:
            ax.axvline(p, color=INK, lw=1.2, ymin=1.05, ymax=1.25, clip_on=False)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlim(1, sub["position"].max())
        ax.set_ylabel(variant, color=INK, fontsize=10, rotation=0, ha="right", va="center")
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=8)
    axes[-1].set_xlabel("residue position", color=MUTED, fontsize=9)
    axes[0].set_title(f"{PROTEIN_NAME[protein_key]} — dominant (rank-1) SAE feature category "
                       f"along the sequence  (| marks ROI)", color=INK, fontsize=12, loc="left")

    cats_present = [c for c in CATEGORY_ORDER if c in set(rank1["category"])]
    handles = [plt.Line2D([], [], marker="s", ls="", mfc=CATEGORY_COLOR[c], mec="none", ms=10, label=c)
               for c in cats_present]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 0.95),
               frameon=False, fontsize=8, title="category")
    fig.tight_layout()
    out = FIGURES / f"feature_track_hp3_{protein_key}.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protein", choices=["gp8", "tcf"], required=True)
    args = ap.parse_args()
    FIGURES.mkdir(parents=True, exist_ok=True)

    df = scan(args.protein)
    update_catalog(df)
    flag_variant_disagreements(args.protein, df)
    plot_track(args.protein, df)


if __name__ == "__main__":
    main()
