# phage-smallautoencoder

**[Download the results report (.pptx)](https://github.com/espickle1/phage-smallautoencoder/blob/main/manuscript/slides/hp3_tcf_gp8_report.pptx)**

ESMC sparse-autoencoder (SAE) feature analysis of two structural proteins from
*Escherichia* phage HP3 (NC_041919.1): the baseplate wedge subunit **gp8**
(YP_010228805.1 / evolved variant URQ01383.1) and the **tail collar fiber**
(YP_010228801.1 / evolved variant URQ01387.1).

For each protein, sequences are read residue-by-residue through ESMC-6B and a
layer-60, 16,384-feature sparse autoencoder. The repo tracks two labeled
sequence variants per protein — referred to throughout as "wt" and "hp3e" (or
"evolved") — and compares their SAE-feature and AlphaFold-3 structure profiles
at a defined set of residue positions.

Per `data/hp3_to_hp3e.txt` and `data/gp8_alignment.csv`:

- **Tail collar fiber**: wt and hp3e are both 516 aa, differing at three
  positions — Y444H, K464R, Y465H.
- **gp8**: wt is 334 aa, hp3e is 337 aa; the difference is a 3-residue DPN
  tandem duplication at wt positions 287–290.

## Repository layout

- `src/` — analysis and figure-generation scripts.
  - `hp3_feature_cache.py`, `hp3_feature_lookup.py` — SAE feature ID → text
    description lookups.
  - `feature_matrix_hp3.py`, `hp3_whole_protein_features.py`,
    `mutation_heatmap_hp3.py`, `residue_map_hp3.py`, `residue_clusters_hp3.py`
    — per-residue SAE/embedding feature extraction, clustering, and plotting.
  - `cluster_structure_hp3.py`, `roi_structure_figures_hp3.py`,
    `gp8_neighbor_context_hp3.py`, `tcf_interface_contacts.py`,
    `tcf_trimer_figures.py`, `annotate_t4_overview.py` — AlphaFold-3 structure
    rendering and geometric analysis (ChimeraX `.cxc` scripts, contact-distance
    measurements).
  - `make_hp3_tcf_gp8_report.py`, `make_hp3_dpn_slides.py` — build the
    `.pptx` report decks in `manuscript/slides/` from the figures and CSVs
    above.
  - `colab_export_hp3.py` — Colab cell that exports raw ESMC-6B SAE
    activations, layer-60 embeddings, and masked-marginal amino-acid
    probabilities for the four sequences (gp8 wt/hp3e, tail collar fiber
    wt/hp3e) to `analysis/features/` and `analysis/structures/`.
- `data/` — input sequences (`sequences/`), the wt↔hp3e alignments/diffs, and
  AlphaFold-3 structure predictions (`structures/`).
- `analysis/` — derived outputs: per-residue feature tables and clusters
  (`features/`), rendered figures (`figures/`), and structure-derived contact
  data (`structures/`).
- `manuscript/slides/` — generated `.pptx` report decks and the Baylor logo
  asset used in their title slides.

## Reproducing

Each `src/*.py` script resolves paths relative to the repo root
(`Path(__file__).resolve().parent.parent`) and can be run directly, e.g.:

```
python src/make_hp3_tcf_gp8_report.py            # 16:9 deck
ASPECT=4:3 python src/make_hp3_tcf_gp8_report.py # 4:3 deck
python src/make_hp3_dpn_slides.py
```

Dependencies: `esm` (Biohub/ESMC client), `torch`, `numpy`, `pandas`,
`scipy`, `scikit-learn`, `gemmi`, `matplotlib`, `adjustText`, `Pillow`,
`python-pptx`. Structure rendering scripts additionally shell out to
[ChimeraX](https://www.cgl.ucsf.edu/chimerax/) via its `.cxc` scripts.

Raw SAE/embedding features (`analysis/features/*.npz`) and masked-marginal
probabilities (`analysis/structures/masked_marginals_hp3.csv`) are produced
by `src/colab_export_hp3.py`, run in a Colab/Jupyter session with an
authenticated ESMC-6B client — see the docstring at the top of that file for
the exact steps.

## License

Apache License 2.0 — see `LICENSE`.
