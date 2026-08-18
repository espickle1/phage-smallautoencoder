"""Heatmaps of ESMC masked-marginal amino-acid preference at the HP3 gp8 and
tail_collar_fiber ROIs (own script, own data file, own figures).

Reads ``analysis/structures/masked_marginals_hp3.csv`` (from
``colab_export_hp3.py``): one row per ROI position with the 20
amino-acid probabilities obtained by masking that position and predicting from
context.

The two proteins here have genuinely different ROI shapes:
  - tail_collar_fiber: 3 ordinary point substitutions (Y444H, K464R, Y465H),
    masked in WT context.
  - gp8: a 3-residue "DPN" insertion, which has no single WT position to
    substitute. Rows are masked in WT context (the flank just before the
    insertion) or EVOLVED context (the 3 inserted residues themselves, plus
    the flank just after) -- see colab_export_hp3.py's module docstring for
    the full design rationale.

Because every row's ``original`` column holds the amino acid actually present
in whichever sequence supplied the masking context (WT AA for substitutions
and the WT flank; the evolved AA itself for insertion/evolved-flank rows),
the "boxed" cell and the log P(aa)/P(original) LLR are well-defined for every
row, substitution or not. Only substitution rows additionally get a circled
``variant`` marker.

Two views (``--mode``):
  prob : P(amino acid | masked context) - sequential blue.
  llr  : log( P(aa) / P(original) ) - diverging; blue = model favours over
         the AA actually present, red = disfavours.

    python mutation_heatmap_hp3.py --mode both     # default
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent.parent
CSV = ROOT / "analysis" / "structures" / "masked_marginals_hp3.csv"
FIG = ROOT / "analysis" / "figures"

AAS = list("ACDEFGHIKLMNPQRSTVWY")
Y_ORDER = list("KRH" "DE" "NQST" "CGP" "AVILMFYW")  # positive/negative/polar/special/hydrophobic

BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
             "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
PROB_CMAP = LinearSegmentedColormap.from_list("blues_ramp", BLUE_RAMP)
DIV_CMAP = LinearSegmentedColormap.from_list("rd_gy_bu", [
    "#7f1d1d", "#c0392b", "#e34948", "#efb0af", "#f0efec",
    "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
INK, MUTED = "#0b0b0b", "#898781"
PROTEIN_COLOR = {"tcf": "#1baf7a", "gp8": "#2a78d6"}

# Explicit row order (protein, context_variant, position): tcf's 3 substitutions
# first, then gp8's story left-to-right across the insertion junction (WT flank
# -> the 3 inserted residues -> evolved flank).
ORDER = [
    ("tcf", "wt", 444), ("tcf", "wt", 464), ("tcf", "wt", 465),
    ("gp8", "wt", 290),
    ("gp8", "evolved", 290), ("gp8", "evolved", 291), ("gp8", "evolved", 292),
    ("gp8", "evolved", 293),
]

KIND_TAG = {
    ("tcf", "wt"): "",
    ("gp8", "wt"): " wt-flank",
    ("gp8", "evolved"): "",  # disambiguated per-row below (insertion vs evolved-flank)
}


def load() -> pd.DataFrame:
    if not CSV.exists():
        raise FileNotFoundError(
            f"{CSV} not found. Run handoff/colab_export_hp3.py on Colab first, "
            f"then download masked_marginals_hp3.csv into analysis/structures/."
        )
    df = pd.read_csv(CSV)
    order_idx = {key: i for i, key in enumerate(ORDER)}
    df["_o"] = [order_idx[(p, c, pos)] for p, c, pos in
                zip(df["protein"], df["context_variant"], df["position"])]
    df = df.sort_values("_o").reset_index(drop=True)
    df["variant"] = df["variant"].where(df["variant"].notna(), "")
    # gp8 position 293 is the evolved-side flank, not an insertion residue --
    # tag it distinctly from the 290-292 insertion rows for the x-axis label.
    df["kind_label"] = [
        KIND_TAG.get((p, c), "") if not (p == "gp8" and c == "evolved" and pos == 293)
        else " ev-flank"
        for p, c, pos in zip(df["protein"], df["context_variant"], df["position"])
    ]
    df["kind_label"] = [
        lbl if not (p == "gp8" and c == "evolved" and pos in (290, 291, 292)) else " ins"
        for p, c, pos, lbl in zip(df["protein"], df["context_variant"], df["position"],
                                   df["kind_label"])
    ]
    return df


def _draw(df: pd.DataFrame, M: np.ndarray, *, imshow_kw: dict, cbar_label: str,
          title: str, out: Path) -> None:
    """Shared layout: matrix M is [20 rows (Y_ORDER) x n residues]."""
    n = len(df)
    labels = [f"{r.protein} {r.original}{r.position}{r.kind_label}" for r in df.itertuples()]
    ypos = {a: i for i, a in enumerate(Y_ORDER)}

    fig, ax = plt.subplots(figsize=(max(8, 0.6 * n + 3), 8))
    im = ax.imshow(M, aspect="auto", **imshow_kw)
    ax.set_xticks(range(n), labels, rotation=90, fontsize=8)
    ax.set_yticks(range(len(Y_ORDER)), list(Y_ORDER), fontsize=9)
    ax.set_ylabel("amino acid", color=MUTED, fontsize=10)
    for tick, prot in zip(ax.get_xticklabels(), df["protein"]):
        tick.set_color(PROTEIN_COLOR.get(prot, INK))

    for j, r in enumerate(df.itertuples()):
        ax.add_patch(Rectangle((j - 0.5, ypos[r.original] - 0.5), 1, 1,
                               fill=False, edgecolor=INK, lw=1.6, zorder=4))
        if isinstance(r.variant, str) and r.variant in ypos:
            ax.plot(j, ypos[r.variant], marker="o", ms=5, mfc="none",
                    mec="#111111", mew=1.6, zorder=5)

    bounds = df.groupby("protein", sort=False)["position"].count().cumsum().tolist()
    starts = [0] + bounds[:-1]
    for b in bounds[:-1]:
        ax.axvline(b - 0.5, color=INK, lw=1.2, zorder=6)
    for s, e, prot in zip(starts, bounds, df["protein"].unique()):
        ax.text((s + e - 1) / 2, -0.85, prot, ha="center", va="bottom",
                fontsize=11, color=INK, fontweight="bold")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(cbar_label, color=INK, fontsize=9)
    ax.set_title(title, color=INK, fontsize=12, loc="left", pad=28)
    handles = [
        plt.Line2D([], [], marker="s", ls="", mfc="none", mec=INK, mew=1.6, ms=10,
                   label="AA present in this context"),
        plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#111111", mew=1.6, ms=8,
                   label="named variant AA (substitutions only)"),
    ] + [plt.Line2D([], [], marker="s", ls="", color=c, ms=9, label=p)
         for p, c in PROTEIN_COLOR.items()]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.12, 1.0),
              frameon=False, fontsize=8, title="labels")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _matrix(vals: np.ndarray) -> np.ndarray:
    """[n, 20] in AAS order -> [20, n] in Y_ORDER rows."""
    return vals[:, [AAS.index(a) for a in Y_ORDER]].T


def render_prob(df: pd.DataFrame) -> None:
    M = _matrix(df[AAS].to_numpy(dtype=float))
    _draw(df, M, imshow_kw=dict(cmap=PROB_CMAP, vmin=0, vmax=float(M.max())),
          cbar_label="P(amino acid | masked context)",
          title="ESMC masked-marginal amino-acid probability at HP3 gp8 / tail_collar_fiber ROIs",
          out=FIG / "mutation_heatmap_hp3.png")


def render_llr(df: pd.DataFrame) -> None:
    prob = df[AAS].to_numpy(dtype=float)
    p_orig = np.array([prob[i, AAS.index(o)] for i, o in enumerate(df["original"])])
    llr = np.log(prob / p_orig[:, None])          # ln P(aa)/P(original); original cell = 0
    M = _matrix(llr)
    vlim = float(np.quantile(np.abs(llr), 0.95))
    _draw(df, M, imshow_kw=dict(cmap=DIV_CMAP, vmin=-vlim, vmax=vlim),
          cbar_label=f"log P(aa) / P(context AA), clipped +/-{vlim:.1f}   —   "
                     f"blue: favoured over the AA present",
          title="ESMC preference vs the AA present in context (log-likelihood ratio) — HP3",
          out=FIG / "mutation_heatmap_hp3_llr.png")


def summary(df: pd.DataFrame) -> None:
    print(f"{'residue':22} {'AA':>3} {'p(AA)':>6} {'var':>3} {'p(var)':>7} "
          f"{'LLR(var)':>8} {'argmax':>6} {'p(max)':>6}")
    rows = []
    for r in df.itertuples():
        p = {a: getattr(r, a) for a in AAS}
        top = max(p, key=p.get)
        has_var = isinstance(r.variant, str) and r.variant in p
        llr_var = np.log(p[r.variant] / p[r.original]) if has_var else np.nan
        name = f"{r.protein} {r.original}{r.position}{r.kind_label}"
        rows.append((name, r.original, p[r.original],
                     r.variant if has_var else "-",
                     f"{p[r.variant]:.3f}" if has_var else "  -  ",
                     f"{llr_var:+.2f}" if has_var else "  -  ", top, p[top]))
    def _key(row):
        return -float(row[5]) if row[5].strip() != "-" else 1e9
    for name, wt, pwt, var, pvar, llr, top, pmax in sorted(rows, key=_key):
        print(f"{name:22} {wt:>3} {pwt:6.3f} {var:>3} {pvar:>7} {llr:>8} {top:>6} {pmax:6.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mode", choices=["prob", "llr", "both"], default="both")
    args = ap.parse_args()
    df = load()
    if args.mode in ("prob", "both"):
        render_prob(df)
    if args.mode in ("llr", "both"):
        render_llr(df)
    summary(df)


if __name__ == "__main__":
    main()
