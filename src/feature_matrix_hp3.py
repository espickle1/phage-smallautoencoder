"""Feature x ROI matrix for HP3 gp8 / tail_collar_fiber -- the HP3 analog of
``feature_matrix.py``, run one protein at a time (own figure per protein, per
the "fully separate analyses" decision -- no combined gp8+tcf figure).

Reads ``analysis/features/hp3_roi_features_top5.csv`` (from
``src/hp3_feature_lookup.py``) and draws a sparse dot matrix:

* x-axis = that protein's ROI positions, grouped by variant (wt -> evolved).
* y-axis = every distinct feature in the table, most-shared (across ROIs) at top.
* a dot at (ROI, feature) means that feature is in the ROI's top-5; dot size
  scales with activation, dot colour encodes the feature's functional category.

Output:
    analysis/figures/feature_matrix_hp3_<protein>.png

Examples
--------
    python feature_matrix_hp3.py --protein gp8
    python feature_matrix_hp3.py --protein tcf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
FEATURES = ROOT / "analysis" / "features"
FIGURES = ROOT / "analysis" / "figures"

PROTEIN_NAME = {"gp8": "gp8", "tcf": "tail_collar_fiber"}
VARIANT_ORDER = ["wt", "evolved"]
VARIANT_COLOR = {"wt": "#2a78d6", "evolved": "#eda100"}

CATEGORY_ORDER = [
    "Domain", "Disorder", "Interaction site", "Catalytic function",
    "Compositional bias", "Ligand-binding site", "Structural motif", "Uncategorized",
]
CATEGORY_COLOR = {
    "Domain": "#2a78d6", "Disorder": "#1baf7a", "Interaction site": "#eda100",
    "Catalytic function": "#008300", "Compositional bias": "#4a3aa7",
    "Ligand-binding site": "#e34948", "Structural motif": "#e87ba4",
    "Uncategorized": "#898781",
}
INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"

SIZE_FLOOR, SIZE_SCALE = 24.0, 400.0


def _size(activation) -> float:
    return SIZE_FLOOR + SIZE_SCALE * activation


def load(protein_key: str) -> pd.DataFrame:
    path = FEATURES / "hp3_roi_features_top5.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run src/hp3_feature_lookup.py first.")
    df = pd.read_csv(path)
    df = df[df["protein"] == PROTEIN_NAME[protein_key]].reset_index(drop=True)
    df["roi"] = [f"{v} {r}{int(p)}" for v, r, p in zip(df.variant, df.residue, df.position)]
    df["category"] = df["category"].where(df["category"].isin(CATEGORY_COLOR), "Uncategorized")
    return df


def draw(df: pd.DataFrame, out: Path, protein_key: str) -> None:
    rois = (df[["variant", "position", "roi"]].drop_duplicates()
            .assign(_o=lambda d: d.variant.map(VARIANT_ORDER.index))
            .sort_values(["_o", "position"]))
    roi_labels = rois["roi"].tolist()
    roi_x = {r: i for i, r in enumerate(roi_labels)}

    feat = (df.groupby("feature_id")
              .agg(label=("label", "first"), category=("category", "first"),
                   n_rois=("roi", "nunique"), max_act=("activation", "max"))
              .sort_values(["n_rois", "max_act"], ascending=[False, False]))
    feat_labels = [f"{r.label} ({fid})" for fid, r in feat.iterrows()]
    feat_y = {fid: i for i, fid in enumerate(feat.index)}

    n_roi, n_feat = len(roi_labels), len(feat)
    fig, ax = plt.subplots(figsize=(0.6 * n_roi + 7.5, 0.32 * n_feat + 2.5))

    ax.set_axisbelow(True)
    ax.grid(color=GRID, lw=0.5, zorder=0)

    for r in df.itertuples():
        ax.scatter(roi_x[r.roi], feat_y[r.feature_id], s=_size(r.activation),
                   c=CATEGORY_COLOR[r.category], edgecolors="white", linewidths=0.6, zorder=3)

    ax.set_xticks(range(n_roi), roi_labels, rotation=90, fontsize=8, color=INK)
    ax.set_yticks(range(n_feat), feat_labels, fontsize=7.5, color=INK)
    ax.set_xlim(-0.5, n_roi - 0.5)
    ax.set_ylim(-0.5, n_feat - 0.5)
    ax.invert_yaxis()
    ax.tick_params(colors=MUTED, length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    counts = rois.groupby("variant", sort=False)["roi"].count()
    counts = counts.reindex([v for v in VARIANT_ORDER if v in counts.index])
    bounds = counts.cumsum().tolist()
    starts = [0] + bounds[:-1]
    for b in bounds[:-1]:
        ax.axvline(b - 0.5, color=INK, lw=1.0, zorder=1)
    for s, e, variant in zip(starts, bounds, counts.index):
        ax.text((s + e - 1) / 2, -1.0, variant, ha="center", va="bottom",
                fontsize=11, fontweight="bold", color=VARIANT_COLOR[variant])

    cats_present = [c for c in CATEGORY_ORDER if c in set(df["category"])]
    cat_handles = [plt.Line2D([], [], marker="o", ls="", mfc=CATEGORY_COLOR[c],
                              mec="white", mew=0.6, ms=9, label=c) for c in cats_present]
    leg1 = ax.legend(handles=cat_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                     frameon=False, fontsize=8, title="feature category",
                     handletextpad=0.6, labelspacing=0.7)
    leg1.get_title().set_color(INK)
    ax.add_artist(leg1)

    size_vals = [0.25, 0.5, 1.0]
    size_handles = [plt.Line2D([], [], marker="o", ls="", mfc=MUTED, mec="white",
                               mew=0.6, ms=(_size(v) ** 0.5), label=f"{v:.2f}")
                    for v in size_vals]
    leg2 = ax.legend(handles=size_handles, loc="upper left", bbox_to_anchor=(1.02, 0.55),
                     frameon=False, fontsize=8, title="activation", labelspacing=1.3,
                     handletextpad=0.8, borderpad=1.0)
    leg2.get_title().set_color(INK)

    ax.set_title(f"Top-5 SAE features per {PROTEIN_NAME[protein_key]} ROI  —  shared features at top",
                 color=INK, fontsize=12, loc="left", pad=42)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight",
                bbox_extra_artists=(leg1, leg2))
    plt.close(fig)
    print(f"wrote {out}  ({n_feat} features x {n_roi} ROIs)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protein", choices=["gp8", "tcf"], required=True)
    args = ap.parse_args()
    df = load(args.protein)
    draw(df, FIGURES / f"feature_matrix_hp3_{args.protein}.png", args.protein)


if __name__ == "__main__":
    main()
