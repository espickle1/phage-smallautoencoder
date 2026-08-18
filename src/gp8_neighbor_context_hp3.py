#!/usr/bin/env python3
"""SAE-cluster coloring for gp8, in the context of its real structural neighbors --
the T4 bacteriophage cryo-EM baseplate wedge (PDB 9MKB, "Structure of the
bacteriophage T4 portal-neck-tail complex"), a very close homolog of our HP3 gp8.

Inputs (analysis/structures/gp8_overlay/, provided by the user, already
structurally aligned into one coordinate frame -- no further superposition
needed here):
  7c0bb7382004a0cb.A.cif  -- our AF3 gp8 WT monomer (same model as
                             data/structures/fold_yp_010228805_..._wt_model_0.cif,
                             just repositioned to align onto 9MKB)
  1be381ec3af82a1.A.cif   -- our AF3 gp8 evolved (HP3e) monomer, same deal
  9MKB.KN.cif              -- full T4 virion cryo-EM structure (541 chains, 141 MB)

9MKB chain ZE overlays our WT gp8 model almost exactly (CA centroid distance
0.30 A) -- confirming gp8 is the T4 gp8 baseplate-wedge-protein homolog. A
CA-CA contact search (<10 A) from ZE against all 541 other chains finds
exactly 6 true neighbors, matching known T4 wedge stoichiometry (2x gp8 +
gp7 + gp6 per wedge, entity table confirms):
  ZD           -- second gp8 copy (328 aa)
  YC, Yc       -- gp7 baseplate wedge protein (1004 aa each)
  eB, eU, f9   -- gp6 baseplate wedge protein (648 aa each)
All other 534 chains (tail tube/sheath, capsid, portal, fibers, etc.) are
irrelevant to gp8's immediate wedge context and are dropped entirely (not
just hidden) -- per-user decision to defer gp8's place in the whole tail.

Pipeline:
  1. Extract just {ZE, ZD, YC, Yc, eB, eU, f9} from the 141 MB 9MKB.KN.cif into
     a ~3 MB neighbor-context cif (coordinates untouched, still in the aligned
     frame) -- avoids ever loading the full virion in ChimeraX.
  2. Bake SAE-cluster IDs (reusing cluster_structure_hp3.py's cluster map) into
     the B-factor column of the two ALIGNED AF3 gp8 copies (not the originals
     in data/structures/ -- those aren't in this coordinate frame).
  3. Render two figures (WT-in-context, evolved-in-context): gp8 opaque and
     cluster-colored with ROI spheres (as in cluster_structure_hp3.py), the 6
     neighbor chains desaturated/transparent grey-family colors by protein
     identity (gp8-copy / gp7 / gp6), so the neighbors read as context, not focus.

Usage:
    python src/gp8_neighbor_context_hp3.py                 # both variants
    python src/gp8_neighbor_context_hp3.py --variant wt
    python src/gp8_neighbor_context_hp3.py --no-render
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cluster_structure_hp3 import (  # noqa: E402
    CLUSTER_PALETTE, NOISE_COLOR, ROI_HIGHLIGHT_COLOR, REPO,
    _find_chimerax, bake_bfactors, cluster_labels, load_cluster_map, roi_targets,
)

OVERLAY_DIR = REPO / "analysis" / "structures" / "gp8_overlay"
NEIGHBORS_CIF = OVERLAY_DIR / "9mkb_gp8_neighbors.cif"
SOURCE_9MKB = OVERLAY_DIR / "9MKB.KN.cif"
OUT_FIG_DIR = REPO / "analysis" / "figures" / "cluster_structures_hp3"

AF3_ALIGNED_CIF = {"wt": OVERLAY_DIR / "7c0bb7382004a0cb.A.cif",
                    "evolved": OVERLAY_DIR / "1be381ec3af82a1.A.cif"}

# chain name -> (protein identity, display color, legend text). Colors are
# desaturated/greyed relative to CLUSTER_PALETTE so they read as background.
NEIGHBOR_CHAINS = {
    "ZD": ("gp8 (2nd copy)", "#9fb8d9"),
    "YC": ("gp7", "#d9b8cf"),
    "Yc": ("gp7", "#d9b8cf"),
    "eB": ("gp6", "#d9cba8"),
    "eU": ("gp6", "#d9cba8"),
    "f9": ("gp6", "#d9cba8"),
}
NEIGHBOR_TRANSPARENCY = 65  # percent

# Darker legend-swatch colors for the same three neighbor identities -- the
# cartoon colors above are deliberately desaturated/pale for the 3D render and
# are unreadable as small legend text on white.
NEIGHBOR_LEGEND_COLOR = {
    "gp8 (2nd copy)": "#3d6aa8",
    "gp7": "#a8578f",
    "gp6": "#a8873f",
}


def extract_neighbors_cif() -> Path:
    if NEIGHBORS_CIF.exists() and NEIGHBORS_CIF.stat().st_mtime > SOURCE_9MKB.stat().st_mtime:
        print(f"  {NEIGHBORS_CIF.relative_to(REPO)} already up to date, skipping extraction")
        return NEIGHBORS_CIF
    print(f"  extracting {len(NEIGHBOR_CHAINS) + 1} chains from "
          f"{SOURCE_9MKB.name} ({SOURCE_9MKB.stat().st_size / 1e6:.0f} MB)...")
    st = gemmi.read_structure(str(SOURCE_9MKB))
    keep = {"ZE", *NEIGHBOR_CHAINS.keys()}
    new_st = gemmi.Structure()
    new_st.name = st.name
    new_model = gemmi.Model("1")
    for chain in st[0]:
        if chain.name in keep:
            new_model.add_chain(chain.clone())
    new_st.add_model(new_model)
    new_st.setup_entities()
    doc = new_st.make_mmcif_document()
    doc.write_file(str(NEIGHBORS_CIF))
    print(f"  wrote {NEIGHBORS_CIF.relative_to(REPO)} "
          f"({NEIGHBORS_CIF.stat().st_size / 1e6:.1f} MB)")
    return NEIGHBORS_CIF


def build_cxc(variant: str, gp8_cif: Path, cluster_ids: list[int],
              labels: dict[int, str]) -> tuple[Path, Path]:
    def color_for(c: int) -> str:
        return NOISE_COLOR if c == -1 else CLUSTER_PALETTE[c % len(CLUSTER_PALETTE)]

    palette = ":".join(f"{c},{color_for(c)}" for c in sorted(cluster_ids))
    roi = roi_targets("gp8", variant=variant)

    lines = [
        f"# Auto-generated by src/gp8_neighbor_context_hp3.py -- gp8 ({variant}) "
        "in its real T4-wedge structural context (9MKB), SAE-cluster colored.",
        "set bgColor white", "lighting soft", "graphics silhouettes true",
        f"open {NEIGHBORS_CIF}",   # model #1: 6 real neighbor chains
        "hide #1 atoms", "cartoon #1",
    ]
    for chain, (_, color) in NEIGHBOR_CHAINS.items():
        lines.append(f"color #1/{chain} {color}")
    lines.append(f"transparency #1 {NEIGHBOR_TRANSPARENCY} cartoons")

    lines += [
        f"open {gp8_cif}",       # model #2: our cluster-colored, aligned AF3 gp8
        "hide #2 atoms", "cartoon #2",
        f"color byattribute bfactor #2 palette {palette} noValueColor {NOISE_COLOR}",
    ]
    for pos, _ in roi:
        sel = f"#2:{pos}"
        lines += [
            f"show {sel} atoms", f"style {sel} sphere",
            f"color {sel} {ROI_HIGHLIGHT_COLOR}", f"size {sel} atomRadius 1.1",
        ]
    lines += ["windowsize 800 800", "view", "zoom 0.75"]

    lines.append(f'2dlabels create legend_title text "gp8 ({variant}) -- SAE clusters '
                 f'in T4 wedge context" xpos 0.02 ypos 0.94 size 14 color black bold true')

    # --- horizontal legend below the figure, minimum text (matches the
    # cluster_structure_hp3.py ROI-clusters figure style) ------------------
    # Non-noise clusters only; cluster numbering and the "noise" entry are
    # intentionally dropped. Two columns, laid out bottom-up, half the old
    # row spacing.
    legend_items = [(labels.get(c, ""), color_for(c)) for c in sorted(cluster_ids) if c != -1]
    positions = sorted(p for p, _ in roi)
    if positions:
        lo, hi = positions[0], positions[-1]
        contiguous = positions == list(range(lo, hi + 1))
        if lo == hi:
            res_range = str(lo)
        else:
            res_range = f"{lo}-{hi}" if contiguous else ",".join(str(p) for p in positions)
        legend_items.append((f"Residues of Interest: {res_range}", ROI_HIGHLIGHT_COLOR))
    seen_desc = set()
    for chain, (desc, color) in NEIGHBOR_CHAINS.items():
        if desc in seen_desc:
            continue
        seen_desc.add(desc)
        legend_items.append((f"neighbor: {desc}", NEIGHBOR_LEGEND_COLOR.get(desc, color)))

    ncols, col_x, row_step, bottom_y = 2, [0.03, 0.46], 0.0225, 0.03
    nrows = (len(legend_items) + ncols - 1) // ncols
    for idx, (text, color) in enumerate(legend_items):
        row, col = divmod(idx, ncols)
        ypos = round(bottom_y + (nrows - 1 - row) * row_step, 4)
        lines.append(
            f'2dlabels create legend_{idx} text "{text}" '
            f'xpos {col_x[col]} ypos {ypos} size 11 color {color}'
        )

    out_png = OUT_FIG_DIR / f"gp8_{variant}_in_wedge_context.png"
    lines += [f"save {out_png} width 800 height 800 supersample 3", ""]

    OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
    cxc_path = OUT_FIG_DIR / f"render_gp8_{variant}_wedge_context.cxc"
    cxc_path.write_text("\n".join(lines))
    return cxc_path, out_png


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--variant", choices=["wt", "evolved"], default=None,
                     help="default: both")
    ap.add_argument("--chimerax", default=None)
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()

    chimerax = None
    if not args.no_render:
        chimerax = args.chimerax or _find_chimerax()
        if not chimerax:
            sys.exit("ChimeraX not found; pass --chimerax PATH or use --no-render")

    extract_neighbors_cif()
    labels = cluster_labels("gp8", {})

    variants = [args.variant] if args.variant else ["wt", "evolved"]
    for variant in variants:
        print(f"[gp8 {variant}]")
        cluster_map = load_cluster_map("gp8", variant=variant)
        gp8_cif, _chains = bake_bfactors(
            "gp8", cluster_map, AF3_ALIGNED_CIF[variant], f"gp8_{variant}_aligned_overlay")
        cluster_ids = sorted(set(cluster_map.values()))
        cxc_path, out_png = build_cxc(variant, gp8_cif, cluster_ids, labels)
        print(f"  wrote {cxc_path.relative_to(REPO)}")
        if args.no_render:
            continue
        print("  rendering (ChimeraX window will appear briefly)...")
        subprocess.run([chimerax, str(cxc_path)], check=True)
        print(f"  wrote {out_png.relative_to(REPO)}")


if __name__ == "__main__":
    main()
