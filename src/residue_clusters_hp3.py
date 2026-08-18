"""HDBSCAN clustering of HP3 residues in SAE-feature (cosine) space -- run one
protein at a time, per the 'fully separate analyses' decision.

WHY NOT cluster the 2D PCA scatter directly: both proteins' PCA plots show a
dense blob connected to a long, thin, curving arc -- a classic artifact of a
strong sequential/positional gradient dominating the leading PCs, not
necessarily real cluster structure. K-means on just (pca_x, pca_y) would
likely chop that arc into arbitrary segments. Instead:

1. Cluster on the FULL cosine geometry (pairwise cosine distance over the
   16,384-dim SAE activation vectors, L2-normalized rows) via HDBSCAN with a
   precomputed distance matrix -- density-based, no fixed cluster count
   needed, tolerates the arc shape, and flags low-density residues as noise
   (cluster -1) instead of forcing every residue into a group.
2. Interpret each cluster by its shared top SAE features: for every member
   residue, take its top-5 active features (same recipe as
   hp3_feature_lookup.py); the cluster's "identity" is whichever features
   recur across the most member residues, described via the same public
   biohub API (no token needed).
3. Visualize by re-plotting the existing PCA/t-SNE layout, colored by cluster
   instead of variant, plus a print summary of whether WT and evolved land in
   the same cluster at each aligned ROI position.

Usage:
    python src/residue_clusters_hp3.py --protein gp8
    python src/residue_clusters_hp3.py --protein tcf --min-cluster-size 20
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import cosine_distances
from sklearn.preprocessing import normalize

from hp3_feature_cache import get_feature_info

ROOT = Path(__file__).resolve().parent.parent
FEATURES = ROOT / "analysis" / "features"
FIGURES = ROOT / "analysis" / "figures"
STRUCTURES = ROOT / "analysis" / "structures"
SEQ_DIR = ROOT / "data" / "sequences"

N_FEATURES = 16384
VARIANTS = ["wt", "evolved"]
SEQ_FILES = {"gp8": SEQ_DIR / "baseplate_gp8_pair.fasta",
             "tcf": SEQ_DIR / "tail_collar_fiber.fasta"}
PROTEIN_NAME = {"gp8": "gp8", "tcf": "tail_collar_fiber"}

INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"
NOISE_COLOR = "#c9c7c1"
# dataviz categorical palette (light mode), cycled if more clusters appear.
CLUSTER_PALETTE = ["#2a78d6", "#eda100", "#1baf7a", "#e34948", "#4a3aa7",
                   "#e87ba4", "#008300", "#8b5a2b", "#00a3a3", "#a35200"]


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


def build_meta(protein_key: str) -> pd.DataFrame:
    records = _parse_fasta(SEQ_FILES[protein_key])
    seqs = {"wt": records[0][1], "evolved": records[1][1]}
    roi_keys = _roi_keys(protein_key)
    frames = []
    for variant in VARIANTS:
        seq = seqs[variant]
        df = pd.DataFrame({"position": np.arange(1, len(seq) + 1), "residue": list(seq)})
        df.insert(0, "protein", protein_key)
        df.insert(1, "variant", variant)
        frames.append(df)
    meta = pd.concat(frames, ignore_index=True)
    meta["is_roi"] = [(v, p) in roi_keys for v, p in zip(meta["variant"], meta["position"])]
    meta["label"] = [
        f"{v} {r}{int(p)}" if roi else ""
        for v, r, p, roi in zip(meta["variant"], meta["residue"], meta["position"], meta["is_roi"])
    ]
    return meta


def load_sae(protein_key: str, meta: pd.DataFrame) -> np.ndarray:
    X = np.zeros((len(meta), N_FEATURES), dtype=np.float32)
    off = 0
    for variant in VARIANTS:
        f = np.load(FEATURES / f"{protein_key}_{variant}.npz")
        n = int(f["seqlen"])
        X[off:off + n] = f["features"]
        off += n
    return X


def cluster(protein_key: str, min_cluster_size: int) -> None:
    meta = build_meta(protein_key)
    X = load_sae(protein_key, meta)
    Xn = normalize(X, norm="l2", axis=1)

    print(f"[{protein_key}] {len(meta)} residues; computing {len(meta)}x{len(meta)} "
          f"cosine distance matrix...")
    D = cosine_distances(Xn).astype(np.float64)
    np.fill_diagonal(D, 0.0)

    hdb = HDBSCAN(min_cluster_size=min_cluster_size, metric="precomputed")
    labels = hdb.fit_predict(D)
    meta["cluster"] = labels

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int((labels == -1).sum())
    print(f"  {n_clusters} clusters found; {n_noise}/{len(meta)} residues unclustered (noise)")
    for c in sorted(set(labels)):
        if c == -1:
            continue
        print(f"  cluster {c}: {int((labels == c).sum())} residues")

    out_meta = FEATURES / f"residue_clusters_hp3_{protein_key}.csv"
    meta.to_csv(out_meta, index=False)
    print(f"  wrote {out_meta.relative_to(ROOT)}")

    # --- per-cluster marker features: detection rate + specificity vs. rest ---
    # Single-cell-style marker detection (SAE feature activation per residue is a
    # close analog of gene expression per cell): a feature only counts as
    # "representative" of a cluster if it's (a) consistently active across most
    # of the cluster's members, AND (b) enriched there relative to the rest of
    # the protein -- plain top-activation-per-residue can surface features that
    # are just generically common everywhere, which isn't distinctive.
    MIN_DETECTION_RATE = 0.5
    EPS = 1e-4
    rows = []
    for c in sorted(set(labels)):
        if c == -1:
            continue
        member_idx = np.where(labels == c)[0]
        other_idx = np.where(labels != c)[0]
        mean_in = X[member_idx].mean(axis=0)
        mean_out = X[other_idx].mean(axis=0)
        detection_rate = (X[member_idx] > 0).mean(axis=0)
        specificity = np.log2((mean_in + EPS) / (mean_out + EPS))

        candidates = np.where(detection_rate >= MIN_DETECTION_RATE)[0]
        candidates = candidates[np.argsort(specificity[candidates])[::-1]][:5]

        for rank, fid in enumerate(candidates, 1):
            fid = int(fid)
            info = get_feature_info(fid)
            rows.append(dict(
                cluster=c, n_members=len(member_idx), rank=rank, feature_id=fid,
                detection_rate=round(float(detection_rate[fid]), 3),
                mean_in=round(float(mean_in[fid]), 4), mean_out=round(float(mean_out[fid]), 4),
                specificity_log2fc=round(float(specificity[fid]), 3),
                label=info.get("label", ""), category=info.get("category", "Uncategorized"),
            ))
    feat_df = pd.DataFrame(rows)
    out_feat = FEATURES / f"cluster_features_hp3_{protein_key}.csv"
    feat_df.to_csv(out_feat, index=False)
    print(f"  wrote {out_feat.relative_to(ROOT)}")
    if not feat_df.empty:
        print("\n  cluster identities (top marker feature):")
        for c in sorted(feat_df["cluster"].unique()):
            top = feat_df[(feat_df.cluster == c) & (feat_df["rank"] == 1)].iloc[0]
            print(f"    cluster {c} (n={top.n_members}): {top.label!r} "
                  f"(detected in {top.detection_rate:.0%} of members, "
                  f"log2FC={top.specificity_log2fc:+.2f} vs rest, {top.category})")

    # --- ROI cluster-membership check: does WT vs evolved land in the same cluster? ---
    roi = meta[meta["is_roi"]].sort_values(["position", "variant"])
    if not roi.empty:
        print("\n  ROI cluster membership:")
        for pos, grp in roi.groupby("position"):
            tags = ", ".join(f"{r.variant}:{r.residue}{pos}->cluster {r.cluster}"
                              for r in grp.itertuples())
            same = grp["cluster"].nunique() == 1
            flag = "" if same else "  <-- DIFFERS across variant"
            print(f"    pos {pos}: {tags}{flag}")

    # --- visualize: re-layout in 2D, colored by cluster ---
    pca = PCA(n_components=2, random_state=0).fit_transform(Xn)
    svd50 = TruncatedSVD(n_components=min(50, X.shape[1] - 1), random_state=0).fit_transform(Xn)
    tsne = TSNE(n_components=2, metric="cosine", init="pca", perplexity=30,
                random_state=0).fit_transform(svd50)

    for coords, tag in ((pca, "pca"), (tsne, "tsne")):
        fig, ax = plt.subplots(figsize=(9, 8))
        for c in sorted(set(labels)):
            m = labels == c
            color = NOISE_COLOR if c == -1 else CLUSTER_PALETTE[c % len(CLUSTER_PALETTE)]
            alpha = 0.35 if c == -1 else 0.7
            lbl = "noise" if c == -1 else f"cluster {c}"
            ax.scatter(coords[m, 0], coords[m, 1], s=10, c=color, alpha=alpha,
                       linewidths=0, label=lbl)
        roi_mask = meta["is_roi"].to_numpy()
        ax.scatter(coords[roi_mask, 0], coords[roi_mask, 1], s=180, marker="*",
                   facecolors="none", edgecolors=INK, linewidths=1.3, zorder=5, label="ROI")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_title(f"{PROTEIN_NAME[protein_key]} residues — SAE space ({tag}), "
                     f"HDBSCAN clusters", color=INK, fontsize=12, loc="left")
        ax.set_xlabel("dim 1", color=MUTED, fontsize=9)
        ax.set_ylabel("dim 2", color=MUTED, fontsize=9)
        ax.legend(frameon=False, fontsize=8, loc="best", ncol=2)
        fig.tight_layout()
        out = FIGURES / f"residue_clusters_hp3_{protein_key}_{tag}.png"
        fig.savefig(out, dpi=200, facecolor="white")
        plt.close(fig)
        print(f"  wrote {out.relative_to(ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protein", choices=["gp8", "tcf"], required=True)
    ap.add_argument("--min-cluster-size", type=int, default=15)
    args = ap.parse_args()
    FIGURES.mkdir(parents=True, exist_ok=True)
    cluster(args.protein, args.min_cluster_size)


if __name__ == "__main__":
    main()
