#!/usr/bin/env python3
"""Render ChimeraX structure figures for the HP3 gp8 / tail_collar_fiber ROIs
(own script, own data, own output folder).

Two figure types, both driven by ``analysis/structures/masked_marginals_hp3.csv``
(the same ROI source of truth used by ``residue_map_hp3.py`` / ``mutation_heatmap_hp3.py``)
plus ``data/structures/hp3_structure_map.json`` (protein/variant -> AF3 .cif):

1. Per-residue close-ups (one PNG per ROI residue): target side chain red
   ball-and-stick, 6 A neighbors thin grey sticks, on a
   neutral-grey transparent cartoon. tail_collar_fiber's 3 substitutions are
   rendered on BOTH the WT and evolved structure (direct before/after); gp8's
   ROI positions are rendered only where they're actually defined (WT flank at
   290; evolved insertion 290-292 + evolved flank 293) -- there's no WT
   counterpart to render at 291/292 since those residues don't exist in WT.
   -> analysis/figures/roi_structures_hp3/<protein>/<variant>/<label>.png

2. Per-protein WT+evolved overlay: both monomers opened together, superimposed
   by ChimeraX's matchmaker (structural alignment, not sequence identity --
   holds up fine across gp8's 3-residue insertion), semi-transparent cartoons
   colored by variant (WT blue / evolved orange, matching residue_map_hp3.py's
   palette), ROI residues shown as colored spheres in their own variant's
   color.
   -> analysis/figures/roi_structures_hp3/<protein>_overlay.png

CAVEAT (read before trusting the tcf evolved-444 panel/sphere): the evolved
tail_collar_fiber .cif was predicted before the Y444H correction to the FASTA,
so it carries Tyr (Y) at 444, not the confirmed His (H) -- see the "caveat"
field in hp3_structure_map.json. That panel and sphere are flagged with a "*"
label and rendered anyway (proceed-with-caveat, not blocked on a re-fold).

Usage:
    python src/roi_structure_figures_hp3.py                 # generate + render both figure types
    python src/roi_structure_figures_hp3.py --no-render      # only write the .cxc + manifest
    python src/roi_structure_figures_hp3.py --chimerax /path/to/ChimeraX
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
STRUCT_DIR = REPO / "data" / "structures"
STRUCT_MAP_PATH = STRUCT_DIR / "hp3_structure_map.json"
MASKED_MARGINALS = REPO / "analysis" / "structures" / "masked_marginals_hp3.csv"
OUT_DIR = REPO / "analysis" / "figures" / "roi_structures_hp3"

DEFAULT_CHIMERAX_CANDIDATES = [
    r"C:\Program Files\ChimeraX 1.10\bin\ChimeraX.exe",
    r"C:\Program Files\ChimeraX 1.9\bin\ChimeraX.exe",
    "/Applications/ChimeraX-1.11.1.app/Contents/MacOS/ChimeraX",
]

VARIANT_COLOR = {"wt": "#2a78d6", "evolved": "#eda100"}
CAVEAT_KEYS = {("tail_collar_fiber", "evolved", 444)}


def _find_chimerax() -> str | None:
    for c in DEFAULT_CHIMERAX_CANDIDATES:
        if Path(c).exists():
            return c
    return None


def load_structure_map() -> dict:
    return json.loads(STRUCT_MAP_PATH.read_text())


def load_roi_targets(smap: dict) -> list[dict]:
    """One dict per (protein, variant, position) close-up to render.

    tail_collar_fiber: the 3 substitution positions, on BOTH variants.
    gp8: exactly the (context_variant, position) rows already scored in the
    masked-marginal CSV (no synthetic "same position, other variant" rows --
    positions past the insertion junction don't have a real counterpart).

    ``residue`` is always read from that variant's OWN structure (CA atoms),
    not copied from the masked-marginals CSV -- the CSV's ``original`` column
    is the AA in whichever sequence supplied the *masking context* (always
    WT for tcf's substitutions), so blindly reusing it for the evolved variant
    would mislabel evolved residues with their WT identity (e.g. tag evolved
    464 "K464" instead of "R464"). Reading straight from each structure's CA
    atoms is also how the tcf-evolved-444 caveat (Y in this structure, not the
    confirmed H) surfaces automatically, without special-casing it here.
    """
    mm = pd.read_csv(MASKED_MARGINALS)
    ca_cache: dict[str, dict[int, str]] = {}

    def residue_at(protein_key: str, variant: str, pos: int) -> str:
        key = f"{protein_key}_{variant}"
        if key not in ca_cache:
            ca_cache[key] = parse_cif_ca(STRUCT_DIR / smap[key]["cif"])
        return ca_cache[key][pos]

    targets = []
    for protein_key, group in mm.groupby("protein"):
        protein = {"gp8": "gp8", "tcf": "tail_collar_fiber"}[protein_key]
        if protein_key == "tcf":
            positions = sorted(set(group["position"]))
            for variant in ("wt", "evolved"):
                for pos in positions:
                    targets.append(dict(protein=protein, protein_key=protein_key,
                                         variant=variant, position=int(pos),
                                         residue=residue_at(protein_key, variant, int(pos))))
        else:
            for _, row in group.iterrows():
                variant, pos = row["context_variant"], int(row["position"])
                targets.append(dict(protein=protein, protein_key=protein_key,
                                     variant=variant, position=pos,
                                     residue=residue_at(protein_key, variant, pos)))
    return targets


THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def parse_cif_ca(path: Path) -> dict[int, str]:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        f = line.split()
        if f[3] != "CA":
            continue
        out[int(f[15])] = THREE2ONE.get(f[5], "X")
    return out


def closeup_block(model_n: int, position: int, code: str, out_path: Path) -> str:
    S = f"#{model_n}:{position}"
    env = f"({S} :<6 & #{model_n})"
    return "\n".join([
        f"# --- {code} -> {out_path.name} ---",
        f"hide #{model_n} atoms",
        "~label",
        f"color #{model_n} #bfbfbf",
        f"transparency #{model_n} 55 cartoons",
        f"show {S} atoms",
        f"style {S} ball",
        f"size {S} stickRadius 0.3",
        f"color {S} byhetero",
        f"color {S} & C #d62728",
        f"show {env} & sidechain atoms",
        f"style {env} stick",
        f"size {env} stickRadius 0.12",
        f"color {env} & ~{S} & C #9a9a9a",
        f'label {S} residues text "{code}" height 2.5 bgColor white',
        f"view {env} pad 0.2",
        f"save {out_path} width 1600 height 1200 supersample 3",
        "",
    ])


def build_closeups_cxc(targets: list[dict], smap: dict) -> tuple[str, list[dict]]:
    """One ChimeraX model per (protein,variant), each opened once; close-ups
    for every target on that model reuse it. Returns (cxc_text, manifest_rows)."""
    by_pv = {}
    for t in targets:
        by_pv.setdefault((t["protein_key"], t["variant"]), []).append(t)

    lines = [
        "# Auto-generated by src/roi_structure_figures_hp3.py -- do not edit by hand.",
        "set bgColor white", "lighting soft", "graphics silhouettes true",
    ]
    manifest = []
    model_n = 1
    for (protein_key, variant), items in by_pv.items():
        entry = smap[f"{protein_key}_{variant}"]
        cif_path = STRUCT_DIR / entry["cif"]
        lines.append(f"open {cif_path}")
        lines += [f"hide #{model_n} atoms", f"cartoon #{model_n}", f"color #{model_n} #bfbfbf", ""]
        for t in items:
            pos, protein = t["position"], t["protein"]
            code = f"{t['residue']}{pos}"
            flagged = (protein, variant, pos) in CAVEAT_KEYS
            if flagged:
                code += "*"
            out = OUT_DIR / protein_key / variant / f"{code.rstrip('*')}.png"
            lines.append(closeup_block(model_n, pos, code, out))
            manifest.append(dict(protein=protein, variant=variant, position=pos,
                                  residue=t["residue"], flagged=flagged, png=str(out.relative_to(REPO))))
        model_n += 1
    lines.append("exit")
    return "\n".join(lines), manifest


def build_overlay_cxc(protein_key: str, protein: str, targets: list[dict], smap: dict) -> str:
    wt = smap[f"{protein_key}_wt"]
    ev = smap[f"{protein_key}_evolved"]
    wt_path = STRUCT_DIR / wt["cif"]
    ev_path = STRUCT_DIR / ev["cif"]
    out = OUT_DIR / f"{protein_key}_overlay.png"

    lines = [
        f"# Auto-generated by src/roi_structure_figures_hp3.py -- {protein} WT+evolved overlay.",
        "set bgColor white", "lighting soft", "graphics silhouettes true",
        f"open {wt_path}",   # model #1 = wt
        f"open {ev_path}",   # model #2 = evolved
        "hide #1,2 atoms", "cartoon #1,2",
        "mm #2 to #1",       # structural superposition, not sequence-based
        f"color #1 {VARIANT_COLOR['wt']}",
        f"color #2 {VARIANT_COLOR['evolved']}",
        "transparency #1,2 60 cartoons",
        "",
    ]
    # No per-atom 3D text labels here: the ROI residues sit within a few
    # angstroms of each other (that's the whole point of the figure), so
    # in-scene labels collide illegibly regardless of camera angle. A fixed
    # on-screen legend (below) identifies every sphere unambiguously instead.
    legend = {"wt": [], "evolved": []}
    for t in sorted(targets, key=lambda t: t["position"]):
        variant, pos = t["variant"], t["position"]
        model_n = 1 if variant == "wt" else 2
        color = VARIANT_COLOR[variant]
        code = f"{t['residue']}{pos}"
        if (protein, variant, pos) in CAVEAT_KEYS:
            code += "*"
        S = f"#{model_n}:{pos}"
        lines += [
            f"show {S} atoms", f"style {S} sphere",
            f"color {S} {color}",
        ]
        legend[variant].append(code)
    lines += ["view pad 0.15", ""]  # extra margin so the fold clears the top-left legend

    # ypos empirically verified (2026-07-21): ChimeraX 1.10 2dlabels clip/vanish
    # entirely above ~0.85 and below ~0.10 in this window size -- keep legend
    # text inside that safe band. Top-left corner (not bottom-left): tcf's
    # elongated fiber fold runs bottom-left to top-right and its bottom end
    # crosses right through bottom-left text; top-left is clear in both
    # proteins' renders.
    legend_lines = [
        f"WT ({VARIANT_COLOR['wt']}): {', '.join(legend['wt'])}",
        f"evolved ({VARIANT_COLOR['evolved']}): {', '.join(legend['evolved'])}",
    ]
    for i, (txt, variant) in enumerate(zip(legend_lines, ("wt", "evolved"))):
        lines.append(
            f'2dlabels create legend_{variant} text "{txt}" xpos 0.02 ypos {0.80 - 0.06 * i} '
            f'size 16 color {VARIANT_COLOR[variant]}'
        )
    if any((protein, t["variant"], t["position"]) in CAVEAT_KEYS for t in targets):
        lines.append(
            '2dlabels create legend_caveat '
            'text "* structure predicted before the Y444H fix -- shows Y, not the confirmed H" '
            'xpos 0.02 ypos 0.68 size 14 color black'
        )
    lines += [
        f"save {out} width 1600 height 1200 supersample 3",
        "",
    ]
    return "\n".join(lines), out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chimerax", default=None)
    ap.add_argument("--no-render", action="store_true")
    args = ap.parse_args()

    smap = load_structure_map()
    targets = load_roi_targets(smap)

    for protein_key, entry in smap.items():
        if isinstance(entry, dict) and "cif" in entry:
            cif = STRUCT_DIR / entry["cif"]
            if not cif.exists():
                sys.exit(f"missing structure: {cif}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for pk in ("gp8", "tcf"):
        for v in ("wt", "evolved"):
            (OUT_DIR / pk / v).mkdir(parents=True, exist_ok=True)

    # 1. close-ups
    cxc_close, manifest = build_closeups_cxc(targets, smap)
    cxc_close_path = OUT_DIR / "render_roi_closeups_hp3.cxc"
    cxc_close_path.write_text(cxc_close)
    pd.DataFrame(manifest).to_csv(OUT_DIR / "manifest.csv", index=False)
    print(f"{len(manifest)} close-up panels; wrote {cxc_close_path.relative_to(REPO)}")

    # 2. per-protein overlays
    overlay_paths = []
    for protein_key, protein in [("gp8", "gp8"), ("tcf", "tail_collar_fiber")]:
        prot_targets = [t for t in targets if t["protein_key"] == protein_key]
        cxc_ov, out_png = build_overlay_cxc(protein_key, protein, prot_targets, smap)
        cxc_ov_path = OUT_DIR / f"render_overlay_{protein_key}.cxc"
        cxc_ov_path.write_text(cxc_ov)
        overlay_paths.append((cxc_ov_path, out_png))
        print(f"wrote {cxc_ov_path.relative_to(REPO)} -> {out_png.relative_to(REPO)}")

    if args.no_render:
        print("--no-render: skipping ChimeraX.")
        return

    chimerax = args.chimerax or _find_chimerax()
    if not chimerax or not Path(chimerax).exists():
        sys.exit("ChimeraX not found; pass --chimerax PATH")

    print("Rendering close-ups (a ChimeraX window will appear briefly)...")
    subprocess.run([chimerax, str(cxc_close_path)], check=True)
    for cxc_path, _ in overlay_paths:
        print(f"Rendering {cxc_path.name}...")
        subprocess.run([chimerax, str(cxc_path)], check=True)

    pngs = sorted(OUT_DIR.rglob("*.png"))
    print(f"Done: {len(pngs)} PNGs under {OUT_DIR.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
