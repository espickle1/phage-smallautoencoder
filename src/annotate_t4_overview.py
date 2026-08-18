#!/usr/bin/env python3
"""Annotate the T4 phage schematic for the HP3 report overview slide.

Draws two markers on the public-domain T4 exterior schematic:
  - Tail collar fiber (TCF): a circle on the DISTAL TIP of a projecting collar
    fiber (where the receptor-binding head domain and the ROI 444/464/465 sit),
    NOT on the neck/collar ring where the fiber attaches.
  - Baseplate wedge (gp8): a circle on the baseplate at the bottom of the tail.

Base image: analysis/figures/phage_exterior.png
  (PhageExterior.svg, Wikimedia Commons; author Adenosine; CC BY-SA 3.0;
   rendered to PNG at 1280 px. Attribution required.)

Output: analysis/figures/t4_overview_annotated.png

Usage:
    python src/annotate_t4_overview.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyBboxPatch
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "analysis" / "figures" / "phage_exterior.png"
OUT = ROOT / "analysis" / "figures" / "t4_overview_annotated.png"

TCF_BLUE = "#2A628A"
GP8_BRICK = "#8A3A2A"


def main() -> None:
    base = Image.open(BASE).convert("RGBA")
    W, H = base.size  # 1280 x 1562
    padL, padR, padT, padB = 600, 120, 40, 40
    CW, CH = W + padL + padR, H + padT + padB
    fig = plt.figure(figsize=(CW / 150, CH / 150), dpi=150)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, CW)
    ax.set_ylim(CH, 0)
    ax.axis("off")
    ax.imshow(np.asarray(base), extent=[padL, padL + W, padT + H, padT])

    def img2c(x, y):
        return (x + padL, y + padT)

    FS = 30
    # TCF: circle the distal tip of the left collar/whisker fiber (~486,815 in image px)
    tcf_tip = img2c(486, 815)
    ax.add_patch(Ellipse(tcf_tip, 150, 210, angle=25, fill=False,
                         edgecolor=TCF_BLUE, lw=6))
    tb_x, tb_y, tb_w, tb_h = 60, 470, 500, 150
    ax.add_patch(FancyBboxPatch((tb_x, tb_y), tb_w, tb_h,
                                boxstyle="round,pad=8,rounding_size=18",
                                facecolor=TCF_BLUE, edgecolor="none"))
    ax.text(tb_x + tb_w / 2, tb_y + tb_h / 2, "Tail collar fiber\n(TCF)",
            color="white", fontsize=FS, fontweight="bold", ha="center", va="center")
    ax.plot([tb_x + tb_w, tcf_tip[0] - 80], [tb_y + tb_h / 2, tcf_tip[1]],
            color=TCF_BLUE, lw=5, solid_capstyle="round")

    # gp8: baseplate wedge at the bottom of the tail
    gp8_c = img2c(628, 1315)
    ax.add_patch(Ellipse(gp8_c, 430, 150, fill=False, edgecolor=GP8_BRICK, lw=6))
    gb_w, gb_h = 470, 150
    gb_x, gb_y = CW - gb_w - 30, 1140
    ax.add_patch(FancyBboxPatch((gb_x, gb_y), gb_w, gb_h,
                                boxstyle="round,pad=8,rounding_size=18",
                                facecolor=GP8_BRICK, edgecolor="none"))
    ax.text(gb_x + gb_w / 2, gb_y + gb_h / 2, "Baseplate\nwedge (gp8)",
            color="white", fontsize=FS, fontweight="bold", ha="center", va="center")
    ax.plot([gb_x, gp8_c[0] + 210], [gb_y + gb_h / 2, gp8_c[1]],
            color=GP8_BRICK, lw=5, solid_capstyle="round")

    fig.savefig(OUT, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT.relative_to(ROOT)}  ({Image.open(OUT).size})")


if __name__ == "__main__":
    main()
