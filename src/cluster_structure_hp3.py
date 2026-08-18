#!/usr/bin/env python3
"""Bake HP3 residue-cluster IDs (from residue_clusters_hp3.py) into the B-factor
column of the WT AF3 structure, and render a ChimeraX figure colored by cluster --
so the clusters can be *seen* on the fold instead of read off an ASCII map.

Two outputs per protein:
  analysis/structures/cluster_colored_{protein}_wt.cif  -- same coordinates as the
      source AF3 model, B_iso_or_equiv replaced by cluster ID (-1 = noise/unclustered).
      Portable: any viewer's "color by B-factor" will reproduce the cluster map.
  analysis/figures/cluster_structures_hp3/{protein}_clusters.png -- ChimeraX render,
      cartoon colored with the same discrete palette used by residue_clusters_hp3.py's
      PCA/t-SNE scatters, with a cluster legend (editable via --labels-json once
      the user has renamed clusters).

Uses the WT variant only (both gp8 and tail_collar_fiber structures have contiguous
1..N numbering matching sequence position exactly -- see hp3_structure_map.json).

Usage:
    python src/cluster_structure_hp3.py                       # both proteins
    python src/cluster_structure_hp3.py --protein gp8
    python src/cluster_structure_hp3.py --no-render            # skip ChimeraX
    python src/cluster_structure_hp3.py --labels-json my_labels.json
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
STRUCT_DIR = REPO / "data" / "structures"
STRUCT_MAP_PATH = STRUCT_DIR / "hp3_structure_map.json"
FEATURES = REPO / "analysis" / "features"
OUT_CIF_DIR = REPO / "analysis" / "structures"
OUT_FIG_DIR = REPO / "analysis" / "figures" / "cluster_structures_hp3"

PROTEIN_NAME = {"gp8": "gp8", "tcf": "tail_collar_fiber"}

DEFAULT_CHIMERAX_CANDIDATES = [
    r"C:\Program Files\ChimeraX 1.10\bin\ChimeraX.exe",
    r"C:\Program Files\ChimeraX 1.9\bin\ChimeraX.exe",
    "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX",
]

# Same as residue_clusters_hp3.py, kept in sync so structure colors match the
# PCA/t-SNE scatter colors for the same cluster ID.
NOISE_COLOR = "#c9c7c1"
CLUSTER_PALETTE = ["#2a78d6", "#eda100", "#1baf7a", "#e34948", "#4a3aa7",
                   "#e87ba4", "#008300", "#8b5a2b", "#00a3a3", "#a35200"]


def _find_chimerax() -> str | None:
    for c in DEFAULT_CHIMERAX_CANDIDATES:
        if Path(c).exists():
            return c
    return None


def load_cluster_map(protein_key: str, variant: str = "wt") -> dict[int, int]:
    df = pd.read_csv(FEATURES / f"residue_clusters_hp3_{protein_key}.csv")
    df = df[df["variant"] == variant]
    return dict(zip(df["position"], df["cluster"]))


def resolve_struct_path(protein_key: str, trimer: bool, variant: str,
                         struct_path: str | None) -> Path:
    if struct_path is not None:
        return Path(struct_path)
    smap = json.loads(STRUCT_MAP_PATH.read_text())
    struct_key = f"{protein_key}_3mer_{variant}" if trimer else f"{protein_key}_{variant}"
    return STRUCT_DIR / smap[struct_key]["cif"]


def bake_bfactors(protein_key: str, cluster_map: dict[int, int], src: Path,
                   out_tag: str) -> tuple[Path, list[str]]:
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    missing = 0
    chain_names = []
    # Every chain present is looped over and given the SAME position->cluster
    # map -- correct whenever chains are identical copies (homo-oligomers),
    # which is the only multi-chain case this script has been pointed at so far.
    for chain in st[0]:
        chain_names.append(chain.name)
        for residue in chain:
            cluster = cluster_map.get(residue.seqid.num)
            if cluster is None:
                missing += 1
                cluster = -1
            for atom in residue:
                atom.b_iso = float(cluster)
    if missing:
        print(f"  [{protein_key}] warning: {missing} residue-atoms had no cluster row")

    OUT_CIF_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_CIF_DIR / f"cluster_colored_{out_tag}.cif"
    doc = st.make_mmcif_document()
    doc.write_file(str(out_path))
    print(f"  wrote {out_path.relative_to(REPO)}")
    return out_path, chain_names


ROI_HIGHLIGHT_COLOR = "#000000"

# Hand-tuned camera framing keyed by output tag. Anything not listed keeps the
# default `view` (fit-to-window). Ported from the finalized gp8 2-mer evolved
# ROI figure so a re-render reproduces the framing instead of clobbering it
# (James, 2026-07-23).
ORIENT_OVERRIDES: dict[str, list[str]] = {
    "gp8_2mer_evolved": [
        "cofr #1",
        "zoom 1",
        "move y 5 models #1",
        "turn y 150 models #1",
        "turn x -10 models #1",
    ],
}


def roi_targets(protein_key: str, variant: str = "wt") -> list[tuple[int, str]]:
    """(position, residue) pairs flagged is_roi=True for the given variant."""
    df = pd.read_csv(FEATURES / f"residue_clusters_hp3_{protein_key}.csv")
    df = df[(df["variant"] == variant) & (df["is_roi"])].sort_values("position")
    return list(zip(df["position"].astype(int), df["residue"]))


def cluster_labels(protein_key: str, labels: dict[str, str]) -> dict[int, str]:
    feat = pd.read_csv(FEATURES / f"cluster_features_hp3_{protein_key}.csv")
    top1 = feat[feat["rank"] == 1].set_index("cluster")["label"].to_dict()
    key_prefix = f"{protein_key}:"
    out = {}
    for c, default_label in top1.items():
        out[c] = labels.get(f"{key_prefix}{c}", default_label)
    return out


def build_cxc(protein_key: str, cif_path: Path, cluster_ids: list[int],
              labels: dict[int, str], chain_names: list[str], out_tag: str,
              variant: str = "wt", highlight_roi: bool = False) -> tuple[str, Path]:
    protein = PROTEIN_NAME[protein_key]
    tag = out_tag + ("_roi" if highlight_roi else "")
    out_png = OUT_FIG_DIR / f"{tag}_clusters.png"

    def color_for(c: int) -> str:
        return NOISE_COLOR if c == -1 else CLUSTER_PALETTE[c % len(CLUSTER_PALETTE)]

    palette_terms = ",".join(f"{c}" if False else f"{c},{color_for(c)}"
                              for c in sorted(cluster_ids))
    # ChimeraX byattribute palette: "v1,color1:v2,color2:..." -- since every
    # atom's bfactor is exactly one of these integer control points (no
    # in-between values exist), each atom resolves to an exact color, not an
    # interpolated blend.
    palette = ":".join(f"{c},{color_for(c)}" for c in sorted(cluster_ids))

    label_suffix = f" ({len(chain_names)}-chain)" if len(chain_names) > 1 else ""
    lines = [
        f"# Auto-generated by src/cluster_structure_hp3.py -- {protein}{label_suffix} residues colored by SAE cluster.",
        "set bgColor white", "lighting soft", "graphics silhouettes true",
        f"open {cif_path}",
        "hide #1 atoms", "cartoon #1",
        f"color byattribute bfactor #1 palette {palette} noValueColor {NOISE_COLOR}",
        "view",
        *ORIENT_OVERRIDES.get(out_tag, []),
    ]

    # ROI spheres first, so the ROI legend entry can join the same grid below.
    roi_item = None
    if highlight_roi:
        roi = roi_targets(protein_key, variant=variant)
        positions = sorted(p for p, _ in roi)
        lo, hi = positions[0], positions[-1]
        contiguous = positions == list(range(lo, hi + 1))
        res_range = f"{lo}-{hi}" if contiguous else ",".join(str(p) for p in positions)
        # one selection covering every chain, e.g. #1/A,B:287-293
        sel = f"#1/{','.join(chain_names)}:{res_range}"
        lines += [
            f"show {sel} atoms",
            f"style {sel} sphere",
            f"color {sel} {ROI_HIGHLIGHT_COLOR}",
            f"size {sel} atomRadius 1.1",
        ]
        legend_range = str(lo) if lo == hi else res_range
        roi_item = (f"Residues of Interest: {legend_range}", ROI_HIGHLIGHT_COLOR)

    # --- horizontal legend below the figure ---------------------------------
    # Non-noise clusters only; cluster numbering and the "noise" entry are
    # intentionally dropped (greyed residues read as background on the fold).
    # Two columns, laid out bottom-up so the block sits just under the figure
    # regardless of cluster count. Row spacing is half the old 0.045.
    legend_items = [(labels.get(c, ""), color_for(c))
                    for c in sorted(cluster_ids) if c != -1]
    if roi_item:
        legend_items.append(roi_item)

    ncols, col_x, row_step, bottom_y = 2, [0.03, 0.46], 0.0225, 0.03
    nrows = (len(legend_items) + ncols - 1) // ncols
    last = len(legend_items) - 1
    for idx, (text, color) in enumerate(legend_items):
        row, col = divmod(idx, ncols)
        ypos = round(bottom_y + (nrows - 1 - row) * row_step, 4)
        name = "legend_roi" if roi_item and idx == last else f"legend_{idx}"
        lines.append(
            f'2dlabels create {name} text "{text}" '
            f'xpos {col_x[col]} ypos {ypos} size 11 color {color}'
        )

    lines += [
        f"save {out_png} width 1800 height 1400 supersample 3",
        "",
    ]
    OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
    cxc_path = OUT_FIG_DIR / f"render_{tag}_clusters.cxc"
    cxc_path.write_text("\n".join(lines))
    return str(cxc_path), out_png


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--protein", choices=["gp8", "tcf"], default=None)
    ap.add_argument("--variant", choices=["wt", "evolved"], default="wt",
                     help="which variant's cluster map / ROI flags to apply (default wt)")
    ap.add_argument("--trimer", action="store_true",
                     help="tcf only, and only with the registered structure map: color the "
                          "homotrimer (tcf_3mer_{variant}) instead of the monomer. For any other "
                          "multi-chain structure (e.g. a dimer), use --struct-path instead -- "
                          "every chain found in the file gets the same per-position cluster map.")
    ap.add_argument("--struct-path", default=None,
                     help="path to an arbitrary .cif to color (bypasses hp3_structure_map.json "
                          "lookup entirely). Requires --protein and --tag.")
    ap.add_argument("--tag", default=None,
                     help="output basename tag (used for cluster_colored_{tag}.cif and "
                          "figure/cxc filenames). Required with --struct-path; auto-derived "
                          "otherwise.")
    ap.add_argument("--highlight-roi", action="store_true",
                     help="overlay the observed ROI positions (from residue_clusters_hp3_*.csv's "
                          "is_roi flag) as black spheres, on top of the cluster coloring, "
                          "on every chain")
    ap.add_argument("--chimerax", default=None)
    ap.add_argument("--no-render", action="store_true")
    ap.add_argument("--labels-json", default=None,
                     help='JSON file mapping "protein:cluster" -> renamed label, '
                          'e.g. {"gp8:0": "protease-adjacent linker"}')
    args = ap.parse_args()

    if args.trimer and args.protein != "tcf":
        sys.exit("--trimer only applies to tcf (no gp8 trimer structure)")
    if args.struct_path and not args.protein:
        sys.exit("--struct-path requires --protein (to pick the right cluster map)")
    if args.struct_path and args.trimer:
        sys.exit("--struct-path and --trimer are mutually exclusive -- pass the trimer cif "
                  "directly via --struct-path")

    labels = {}
    if args.labels_json:
        labels = json.loads(Path(args.labels_json).read_text())

    proteins = [args.protein] if args.protein else ["gp8", "tcf"]
    chimerax = None
    if not args.no_render:
        chimerax = args.chimerax or _find_chimerax()
        if not chimerax:
            sys.exit("ChimeraX not found; pass --chimerax PATH or use --no-render")

    for protein_key in proteins:
        trimer = args.trimer and protein_key == "tcf"
        if args.struct_path:
            out_tag = args.tag or Path(args.struct_path).stem
        else:
            out_tag = args.tag or (f"{protein_key}_3mer_{args.variant}" if trimer
                                    else f"{protein_key}_{args.variant}")
        print(f"[{out_tag}]")
        cluster_map = load_cluster_map(protein_key, variant=args.variant)
        src = resolve_struct_path(protein_key, trimer, args.variant, args.struct_path)
        cif_path, chain_names = bake_bfactors(protein_key, cluster_map, src, out_tag)
        cluster_ids = sorted(set(cluster_map.values()))
        labels_for_protein = cluster_labels(protein_key, labels)
        cxc_path, out_png = build_cxc(protein_key, cif_path, cluster_ids, labels_for_protein,
                                       chain_names, out_tag, variant=args.variant,
                                       highlight_roi=args.highlight_roi)
        print(f"  wrote {Path(cxc_path).relative_to(REPO)}")
        if args.no_render:
            continue
        print(f"  rendering (ChimeraX window will appear briefly)...")
        subprocess.run([chimerax, cxc_path], check=True)
        print(f"  wrote {out_png.relative_to(REPO)}")


if __name__ == "__main__":
    main()
