"""Colab cell: export ESMC-6B SAE activations, layer-60 embeddings, and masked-marginal
amino-acid probabilities for the two HP3 WT/evolved protein pairs (gp8, tail_collar_fiber).

WHY THIS SCRIPT IS SELF-CONTAINED: there are two REAL, independently sequenced variants
per protein (WT and evolved), and one of them (gp8) differs by a 3-residue INSERTION, not
a substitution, so there is no single shared "reference + variant list" that a .cif-driven
lookup could use. Sequences are pasted in directly from data/sequences/*.fasta (verified
against the repo copies 2026-07-21), no .cif/get_sequence dependency. AF3 structures for
these four sequences already exist in data/structures/ (see data/structures/hp3_structure_map.json)
for the later ChimeraX ROI-render step -- they are not needed for anything in this script.

HOW TO RUN
1. In a Colab or Jupyter session, set up an authenticated ESMC-6B client with SAE support
   (`ESMCForgeInferenceClient(model="esmc-6b-...", ...)`, requires a Biohub API token) so
   `model` is defined.
2. Paste this whole file into a new cell and run it. Writes, to your working directory (or
   adjust OUT_FEATURES/OUT_STRUCTURES below for a different destination):
     - <name>.npz       dense SAE activations [seqlen, 16384], key "features"
     - <name>_emb.npz   dense layer-60 embeddings [seqlen, 2560], key "emb"
     - masked_marginals_hp3.csv   one row per ROI position x 20 AA columns
   for name in {gp8_wt, gp8_evolved, tcf_wt, tcf_evolved}.
3. Copy those files into the repo:
     <name>.npz, <name>_emb.npz          -> analysis/features/
     masked_marginals_hp3.csv            -> analysis/structures/

ROI DESIGN (read before trusting the masked-marginal rows)
- tail_collar_fiber: three ordinary point substitutions, Y444H, K464R, Y465H (confirmed by direct
  sequence diff of the two 516-aa sequences after a 2026-07-21 manual correction to the FASTA;
  matches data/hp3_to_hp3e.txt -- ignore the FASTA header's vague "99.6% id" claim). Handled by
  masking each position in the WT sequence and reading P(20 AA) from WT context.
- gp8: WT (334 aa) vs evolved (337 aa) differ by a 3-residue "DPN" tandem-repeat INSERTION, not a
  substitution -- confirmed by direct pairwise alignment (see data/gp8_alignment.csv): WT 1-289 ==
  evolved 1-289, evolved gains "DPN" at 290-292, WT 290-334 == evolved 293-337. A substitution-style
  masked-marginal (mask one WT position, compare to one variant AA) doesn't apply to an insertion.
  Instead this script asks two related questions:
    (a) "Does the evolved sequence's own local context consider the inserted residues plausible?"
        -- mask each of evolved positions 290/291/292 (D, P, N) *in evolved context* and read
        P(20 AA) there. High P(D)/P(P)/P(N) would mean the insertion looks unsurprising to the
        model given everything around it; low values would mean the model finds the repeat odd.
    (b) "Does the insertion change what the model expects immediately downstream?" -- mask WT
        position 290 (V, in WT context, no insertion) and the aligned evolved position 293 (V, in
        evolved context, right after the insertion) and compare the two P(20 AA) distributions at
        what is structurally the same residue. A shift here would mean the insertion perturbs local
        context beyond the 3 inserted residues themselves.
  All five gp8 rows are tagged `context` = "wt" or "evolved" in the output so mutation_heatmap.py
  (or its HP3 adaptation) knows which sequence was used, and `kind` = "substitution", "insertion",
  or "flank" so the two protein's rows aren't conflated.
"""

import csv

import numpy as np
import torch
from esm.sdk.api import ESMProtein, LogitsConfig, SAEConfig

try:
    from esm.tokenization import EsmSequenceTokenizer
except Exception:  # pragma: no cover - import path fallback
    from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

# --- output locations (Colab-local by default; point at your Drive folder if you'd rather
# download from there, same pattern as colab_export_embeddings.py) ---
OUT_FEATURES = "."       # <name>.npz, <name>_emb.npz
OUT_STRUCTURES = "."     # masked_marginals_hp3.csv

SAE_LAYER = 60
SAE_MODEL_NAME = "esmc-6b-2024-12-sae-layer60-k64-codebook16384"
AAS = list("ACDEFGHIKLMNPQRSTVWY")

VALID_AAS = set(AAS)


def validate_protein_sequence(sequence: str) -> str:
    sequence = sequence.strip().upper()
    invalid = sorted(set(sequence) - VALID_AAS)
    if invalid:
        raise ValueError(f"Invalid amino acid characters found: {invalid}")
    return sequence


# --- sequences, pasted directly from data/sequences/*.fasta (verified 2026-07-21) ---
SEQUENCES = {
    "gp8_wt": (
        "MNDSSVIYRAIVTSKFRTEKMLNFYNSIGSGPDKNTIFITFGRSEPWSSNENEVGFAPPYPTDSVLGVTD"
        "MWTHMMGTVKVLPSMLDAVIPRRDWGDTRYPDPYTFRINDIVVCNSAPYNATESGAGWLVYRCLDVPDTG"
        "MCSIASLTNKDECLKLGGKWTPSVRSMTPPEGRGDAEGTIEPGDGYVWEYLFEIPPDVSINRCTNEYIVV"
        "PWPEELKEDPTRWGYEDNLTWQQDDFGLIYRVKANTIRFKAYLDSVYFPDAALPGNKGFRQISIITNPLE"
        "AKAHPNDPNVKAEKDYYDPEDLMRHSGEMIYMENRPPIIMAMDQTEEINILFTF"
    ),
    "gp8_evolved": (
        "MNDSSVIYRAIVTSKFRTEKMLNFYNSIGSGPDKNTIFITFGRSEPWSSNENEVGFAPPYPTDSVLGVTD"
        "MWTHMMGTVKVLPSMLDAVIPRRDWGDTRYPDPYTFRINDIVVCNSAPYNATESGAGWLVYRCLDVPDTG"
        "MCSIASLTNKDECLKLGGKWTPSVRSMTPPEGRGDAEGTIEPGDGYVWEYLFEIPPDVSINRCTNEYIVV"
        "PWPEELKEDPTRWGYEDNLTWQQDDFGLIYRVKANTIRFKAYLDSVYFPDAALPGNKGFRQISIITNPLE"
        "AKAHPNDPNDPNVKAEKDYYDPEDLMRHSGEMIYMENRPPIIMAMDQTEEINILFTF"
    ),
    "tcf_wt": (
        "MSNNTYQHVSNESKYVKFDPTGSNFPGTVTTVQSALSKISNIGVNGIPDATMEVKGIAMIASEQEVLDGT"
        "NNSKIVTPATLATRLLYPNATETKYGLTRYSTNEETLKGSDNNSSITPQKLKYHTDDVFKNRYSSESSNG"
        "VIKISSTPAALAGVDDTTAMTPLKTQKLAIKLISQIAPSEDTATESVRGVVQLSTVAQIRQGTLREGYAI"
        "SPYTFMNSVATHEYKGVIRLGTQTEINNNLGGVAVTGETLNGRGATGSMRGVVKLTTQAGIAPEGDSSGA"
        "LAWNADVINTRGGQTINGSLNLDHLTANGIWSRGGMWKNGDQPVATERYASERVPVGTIMMFAGDSAPPG"
        "WIMCHGGTVSGDQFPDYRNVVGTRFGGDWNNPGVPDMRGLFVRGAGTGGHILNQRGQDGYGKDRLGVGCD"
        "GMHVGGVQAQQMSYHKHAGGWGEYNRSEGPFGASVYQGYLGTRKYSDWDNASYFTNDGFELGGPRDAHGT"
        "LNREGLIGYETRPWNISLNYIIKVHY"
    ),
    "tcf_evolved": (
        "MSNNTYQHVSNESKYVKFDPTGSNFPGTVTTVQSALSKISNIGVNGIPDATMEVKGIAMIASEQEVLDGT"
        "NNSKIVTPATLATRLLYPNATETKYGLTRYSTNEETLKGSDNNSSITPQKLKYHTDDVFKNRYSSESSNG"
        "VIKISSTPAALAGVDDTTAMTPLKTQKLAIKLISQIAPSEDTATESVRGVVQLSTVAQIRQGTLREGYAI"
        "SPYTFMNSVATHEYKGVIRLGTQTEINNNLGGVAVTGETLNGRGATGSMRGVVKLTTQAGIAPEGDSSGA"
        "LAWNADVINTRGGQTINGSLNLDHLTANGIWSRGGMWKNGDQPVATERYASERVPVGTIMMFAGDSAPPG"
        "WIMCHGGTVSGDQFPDYRNVVGTRFGGDWNNPGVPDMRGLFVRGAGTGGHILNQRGQDGYGKDRLGVGCD"
        "GMHVGGVQAQQMSYHKHAGGWGEHNRSEGPFGASVYQGYLGTRRHSDWDNASYFTNDGFELGGPRDAHGT"
        "LNREGLIGYETRPWNISLNYIIKVHY"
    ),
}
for _name, _seq in SEQUENCES.items():
    SEQUENCES[_name] = validate_protein_sequence(_seq)
assert len(SEQUENCES["gp8_wt"]) == 334 and len(SEQUENCES["gp8_evolved"]) == 337
assert len(SEQUENCES["tcf_wt"]) == 516 and len(SEQUENCES["tcf_evolved"]) == 516

# ROI targets: (context_sequence_name, position_1indexed, context_aa, variant_aa_or_None, kind)
# context_sequence_name selects which SEQUENCES entry supplies the masking context.
# kind in {substitution, insertion, flank}:
#   substitution -- ordinary point mutation, masked in WT context, compared to the named variant.
#   insertion    -- one of the 3 novel gp8 residues, masked in EVOLVED context (no WT counterpart).
#   flank        -- the residue immediately after the gp8 insertion junction, masked in both WT
#                   and evolved context, to see whether the insertion perturbs local expectation.
TARGETS = [
    ("tcf_wt", 444, "Y", "H", "substitution"),
    ("tcf_wt", 464, "K", "R", "substitution"),
    ("tcf_wt", 465, "Y", "H", "substitution"),
    ("gp8_evolved", 290, "D", None, "insertion"),
    ("gp8_evolved", 291, "P", None, "insertion"),
    ("gp8_evolved", 292, "N", None, "insertion"),
    ("gp8_wt", 290, "V", None, "flank"),
    ("gp8_evolved", 293, "V", None, "flank"),
]


def _to_np(x):
    """ESMC tensors are bfloat16; upcast to float32 on CPU before np.asarray."""
    if isinstance(x, torch.Tensor):
        return x.detach().to(torch.float32).cpu().numpy()
    return np.asarray(x, dtype=np.float32)


def _softmax(x):
    e = np.exp(x - x.max())
    return e / e.sum()


# --- 1. SAE activations + layer-60 embeddings, one pass per sequence ---
for name, seq in SEQUENCES.items():
    protein_tensor = model.encode(ESMProtein(sequence=seq))  # noqa: F821 (notebook global)

    sae_out = model.logits(  # noqa: F821 (notebook global)
        protein_tensor,
        config=LogitsConfig(sae_config=SAEConfig(
            model=SAE_MODEL_NAME, normalize_features=True)),
        return_bytes=False,
    )
    features = sae_out.sae_outputs[SAE_MODEL_NAME].to_dense().numpy()[1:-1]
    assert features.shape[0] == len(seq), (name, features.shape, len(seq))
    np.savez_compressed(f"{OUT_FEATURES}/{name}.npz", features=features, seqlen=len(seq))

    try:
        emb_out = model.logits(  # noqa: F821 (notebook global)
            protein_tensor,
            config=LogitsConfig(return_hidden_states=True, ith_hidden_layer=SAE_LAYER),
            return_bytes=False,
        )
        h = np.squeeze(_to_np(emb_out.hidden_states))
        if h.ndim == 3:
            h = h[SAE_LAYER]
        src = f"hidden_states[layer {SAE_LAYER}]"
    except Exception as e:  # noqa: BLE001 - reason goes in the printout
        emb_out = model.logits(  # noqa: F821 (notebook global)
            protein_tensor, config=LogitsConfig(return_embeddings=True), return_bytes=False)
        h = np.squeeze(_to_np(emb_out.embeddings))
        src = f"embeddings (FINAL layer; layer-{SAE_LAYER} hidden states unavailable: {e})"
    emb = h[1:-1].astype(np.float32)
    assert emb.shape[0] == len(seq), (name, emb.shape, len(seq))
    np.savez_compressed(f"{OUT_FEATURES}/{name}_emb.npz", emb=emb, seqlen=len(seq))
    print(f"{name}: SAE {features.shape}, emb {emb.shape} via {src}")

# --- 2. masked marginals at the ROI positions above ---
tok = EsmSequenceTokenizer()
vocab = tok.get_vocab()
aa_ids = [vocab[a] for a in AAS]
MASK_ID = tok.mask_token_id
print(f"vocab size {len(vocab)}, mask id {MASK_ID}")

rows = []
by_context = {}
for ctx_name, *_ in TARGETS:
    by_context.setdefault(ctx_name, []).append(_)

for ctx_name, targets in by_context.items():
    seq = SEQUENCES[ctx_name]
    protein, variant = ctx_name.split("_")
    pt = model.encode(ESMProtein(sequence=seq))  # noqa: F821 (notebook global)
    base = pt.sequence.clone()
    assert base.shape[0] == len(seq) + 2, (ctx_name, base.shape, len(seq))
    for pos, orig, var, kind in targets:
        assert seq[pos - 1] == orig, (ctx_name, pos, seq[pos - 1], orig)
        pt.sequence[pos] = MASK_ID
        out = model.logits(pt, LogitsConfig(sequence=True), return_bytes=False)  # noqa: F821
        pt.sequence[pos] = base[pos]
        lg = _to_np(out.logits.sequence).squeeze()
        assert max(aa_ids) < lg.shape[-1], (lg.shape, max(aa_ids))
        probs = _softmax(lg[pos][aa_ids])
        rows.append([protein, variant, pos, orig, var or "", kind, *probs.tolist()])
        top3 = sorted(zip(AAS, probs), key=lambda t: -t[1])[:3]
        print(f"  {ctx_name} pos {pos} ({orig}, {kind}) top-3: "
              f"{[(a, round(float(p), 3)) for a, p in top3]}")

out_path = f"{OUT_STRUCTURES}/masked_marginals_hp3.csv"
with open(out_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["protein", "context_variant", "position", "original", "variant",
                "kind", *AAS])
    w.writerows(rows)
print(f"wrote {out_path} ({len(rows)} rows x {len(AAS)} AAs)")
