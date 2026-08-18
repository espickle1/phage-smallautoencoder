"""Residue-similarity maps for the HP3 WT/evolved protein pairs (own script,
own data files, own figures, run one protein at a time -- no pooled
cross-protein basis).

Each *residue* is a point. Distance between residues = how differently they
behave, measured two interchangeable ways (``--rep``):

* ``sae`` — the residue's SAE feature-activation vector (16,384-dim, top-64),
  read locally from ``analysis/features/<protein>_<variant>.npz`` (key
  ``features``).
* ``emb`` — the residue's raw ESMC layer-60 embedding (dense 2,560-dim), read
  from ``analysis/features/<protein>_<variant>_emb.npz`` (key ``emb``).

Cosine metric throughout. For a given ``--protein``, the WT and evolved
sequences are pooled into one shared basis (2 variants), coloured by variant,
with the ROI residues highlighted.
ROI membership is read directly from ``analysis/structures/masked_marginals_hp3.csv``
so it can't drift from the masked-marginal analysis:
  - tail_collar_fiber: the 3 substitution positions (444/464/465), shown in
    BOTH variants (6 ROI points) -- the direct before/after comparison.
  - gp8: the exact (context_variant, position) pairs already scored there --
    WT flank (290), evolved insertion (290-292) + evolved flank (293), 5
    ROI points total. Not every WT position has an evolved counterpart at the
    same index past the insertion, so ROI membership is NOT just "same
    position in both variants" here (see data/gp8_alignment.csv).

Outputs (``--protein gp8 --rep sae`` shown; other combinations swap the names):
    analysis/features/residue_map_coords_hp3_gp8_sae.csv
    analysis/figures/residue_map_hp3_gp8_sae_pca.png
    analysis/figures/residue_map_hp3_gp8_sae_tsne.png

NOTE (2026-07-21): an earlier version also produced an ROI-only pairwise
cosine-distance heatmap + dendrogram (5-6 points). Dropped: with so few ROI
points, each WT/evolved pair is near-trivially its own closest neighbor (same
or near-identical local context), so that figure didn't tell a reader anything
beyond what's already obvious from the ROI definition itself. The full-protein
PCA/t-SNE scatters below are kept because they show the ROI in the context of
every other residue in the protein, which is the actually informative
comparison.

Examples
--------
    python residue_map_hp3.py --protein gp8 --rep sae
    python residue_map_hp3.py --protein tcf --rep emb
    python residue_map_hp3.py --protein gp8 --rep compare
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.manifold import TSNE
from sklearn.preprocessing import normalize

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis"
FEATURES = ANALYSIS / "features"
FIGURES = ANALYSIS / "figures"
STRUCTURES = ANALYSIS / "structures"
SEQ_DIR = ROOT / "data" / "sequences"

N_FEATURES = 16384
VARIANTS = ["wt", "evolved"]
SEQ_FILES = {"gp8": SEQ_DIR / "baseplate_gp8_pair.fasta",
             "tcf": SEQ_DIR / "tail_collar_fiber.fasta"}

# dataviz palette: categorical slots 1-2 (blue / yellow), light mode.
VARIANT_COLOR = {"wt": "#2a78d6", "evolved": "#eda100"}
BLUE_RAMP = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
SEQ_CMAP = LinearSegmentedColormap.from_list("blues_ramp", BLUE_RAMP)
INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """[(header, sequence)] in file order (WT record first, evolved second)."""
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


def _roi_keys(protein: str) -> set[tuple[str, int]]:
    mm = pd.read_csv(STRUCTURES / "masked_marginals_hp3.csv")
    mm = mm[mm["protein"] == protein]
    if protein == "tcf":
        positions = set(mm["position"])
        return {(v, p) for v in VARIANTS for p in positions}
    return set(zip(mm["context_variant"], mm["position"]))


def build_meta(protein: str) -> pd.DataFrame:
    records = _parse_fasta(SEQ_FILES[protein])
    assert len(records) == 2, (protein, len(records))
    seqs = {"wt": records[0][1], "evolved": records[1][1]}
    roi_keys = _roi_keys(protein)

    frames = []
    for variant in VARIANTS:
        seq = seqs[variant]
        df = pd.DataFrame({
            "position": np.arange(1, len(seq) + 1),
            "residue": list(seq),
        })
        df.insert(0, "protein", protein)
        df.insert(1, "variant", variant)
        frames.append(df)
    meta = pd.concat(frames, ignore_index=True)
    meta["is_roi"] = [(v, p) in roi_keys for v, p in zip(meta["variant"], meta["position"])]
    meta["label"] = [
        f"{v} {r}{int(p)}" if roi else ""
        for v, r, p, roi in zip(meta["variant"], meta["residue"], meta["position"], meta["is_roi"])
    ]
    return meta


def load_sae(protein: str, meta: pd.DataFrame) -> np.ndarray:
    X = np.zeros((len(meta), N_FEATURES), dtype=np.float32)
    off = 0
    for variant in VARIANTS:
        path = FEATURES / f"{protein}_{variant}.npz"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not found. Run handoff/colab_export_hp3.py on Colab first "
                f"and place the .npz outputs in analysis/features/."
            )
        f = np.load(path)
        n = int(f["seqlen"])
        X[off:off + n] = f["features"]
        off += n
    return X


def load_emb(protein: str, meta: pd.DataFrame) -> np.ndarray:
    mats = []
    for variant in VARIANTS:
        path = FEATURES / f"{protein}_{variant}_emb.npz"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not found. Run handoff/colab_export_hp3.py on Colab first "
                f"and place the _emb.npz outputs in analysis/features/."
            )
        mats.append(np.load(path)["emb"].astype(np.float32))
    return np.vstack(mats)


REPS = {
    "sae": ("SAE feature activations (16,384-dim, top-64)", load_sae),
    "emb": ("ESMC layer-60 embeddings", load_emb),
}


def _style_axes(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def scatter(protein: str, meta: pd.DataFrame, coords: np.ndarray, title: str, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 8))
    for variant, color in VARIANT_COLOR.items():
        m = (meta["variant"] == variant).to_numpy() & ~meta["is_roi"].to_numpy()
        ax.scatter(coords[m, 0], coords[m, 1], s=8, c=color, alpha=0.45,
                   linewidths=0, label=variant)
    roi = meta["is_roi"].to_numpy()
    for variant, color in VARIANT_COLOR.items():
        m = (meta["variant"] == variant).to_numpy() & roi
        ax.scatter(coords[m, 0], coords[m, 1], s=230, marker="*", c=color,
                   edgecolors=INK, linewidths=1.1, zorder=5)
    texts = [
        ax.text(coords[i, 0], coords[i, 1], meta["label"].iloc[i], fontsize=8,
                color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
        for i in np.where(roi)[0]
    ]
    adjust_text(texts, ax=ax, expand=(1.4, 1.8), force_text=(0.5, 0.8),
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    ax.set_title(title, color=INK, fontsize=12, loc="left")
    ax.set_xlabel("dim 1", color=MUTED, fontsize=9)
    ax.set_ylabel("dim 2", color=MUTED, fontsize=9)
    _style_axes(ax)
    leg = ax.legend(frameon=False, fontsize=9, loc="best", title="variant")
    leg.get_title().set_color(INK)
    fig.tight_layout()
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)


def compare(protein: str, meta: pd.DataFrame) -> None:
    Xs = normalize(load_sae(protein, meta), norm="l2", axis=1)
    Xe = normalize(load_emb(protein, meta), norm="l2", axis=1)

    roi_idx = np.where(meta["is_roi"].to_numpy())[0]
    labels = meta["label"].iloc[roi_idx].to_numpy()
    ds = pdist(Xs[roi_idx], metric="cosine")
    de = pdist(Xe[roi_idx], metric="cosine")
    rho_roi = spearmanr(ds, de).statistic if len(ds) > 1 else float("nan")

    rng = np.random.default_rng(0)
    n = len(meta)
    ii = rng.integers(0, n, 40000)
    jj = rng.integers(0, n, 40000)
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    gs = 1.0 - np.sum(Xs[ii] * Xs[jj], axis=1)
    ge = 1.0 - np.sum(Xe[ii] * Xe[jj], axis=1)
    rho_all = spearmanr(gs, ge).statistic

    fig, ax = plt.subplots(figsize=(6.5, 6))
    variant = meta["variant"].to_numpy()
    pi, pj = np.triu_indices(len(roi_idx), k=1)
    same = variant[roi_idx][pi] == variant[roi_idx][pj]
    ax.scatter(ds[~same], de[~same], s=40, c=MUTED, alpha=0.7, label="wt vs evolved")
    ax.scatter(ds[same], de[same], s=48, c="#2a78d6", alpha=0.85,
               edgecolors=INK, linewidths=0.5, label="same variant")
    ax.plot([0, 1], [0, 1], color=GRID, lw=1, zorder=0)
    ax.set_xlabel("cosine distance — SAE space", color=MUTED, fontsize=9)
    ax.set_ylabel("cosine distance — embedding space", color=MUTED, fontsize=9)
    ax.set_title(f"{protein}: do the two views agree? (ROI pairs)", color=INK,
                 fontsize=12, loc="left")
    ax.text(0.03, 0.97, f"Spearman ρ (ROI pairs) = {rho_roi:.2f}\n"
                        f"Spearman ρ (all residues, sampled) = {rho_all:.2f}",
            transform=ax.transAxes, va="top", fontsize=9, color=INK)
    _style_axes(ax)
    leg = ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    out = FIGURES / f"compare_sae_vs_emb_hp3_{protein}.png"
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"  ROI-pair Spearman rho:      {rho_roi:.3f}")
    print(f"  all-residue Spearman rho:   {rho_all:.3f} (40k sampled pairs)")
    print(f"  wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protein", choices=["gp8", "tcf"], required=True)
    ap.add_argument("--rep", choices=[*REPS, "compare"], default="sae")
    args = ap.parse_args()
    protein, rep = args.protein, args.rep
    FIGURES.mkdir(parents=True, exist_ok=True)

    meta = build_meta(protein)
    print(f"[{protein}] {len(meta)} residues across {len(VARIANTS)} variants; "
          f"{int(meta['is_roi'].sum())} ROI points")

    if rep == "compare":
        print(f"[compare] {protein}: SAE-activation space vs ESMC-embedding space")
        compare(protein, meta)
        return

    desc, loader = REPS[rep]
    print(f"[{rep}] {desc}")
    X = loader(protein, meta)
    print(f"  matrix: {X.shape}")

    Xn = normalize(X, norm="l2", axis=1)
    pca = PCA(n_components=2, random_state=0).fit_transform(Xn)
    svd50 = TruncatedSVD(n_components=min(50, X.shape[1] - 1),
                         random_state=0).fit_transform(Xn)
    tsne = TSNE(n_components=2, metric="cosine", init="pca", perplexity=30,
                random_state=0).fit_transform(svd50)

    scatter(protein, meta, pca, f"{protein} residues — {rep} space (PCA)",
            FIGURES / f"residue_map_hp3_{protein}_{rep}_pca.png")
    scatter(protein, meta, tsne, f"{protein} residues — {rep} space (t-SNE)",
            FIGURES / f"residue_map_hp3_{protein}_{rep}_tsne.png")

    coords = meta.copy()
    coords["pca_x"], coords["pca_y"] = pca[:, 0], pca[:, 1]
    coords["tsne_x"], coords["tsne_y"] = tsne[:, 0], tsne[:, 1]
    coords.to_csv(FEATURES / f"residue_map_coords_hp3_{protein}_{rep}.csv", index=False)
    print(f"  wrote coords + 2 figures for protein={protein} rep={rep}")


if __name__ == "__main__":
    main()
