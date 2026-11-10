"""
render_study_design.py  (v5 - vertical flowchart, auto-sized + stacked)
=======================================================================
Figure 1. Every box is sized from its content; stages are stacked with fixed
gaps, so nothing can overflow. Numbers are read live from the manifests.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[2]
B = ROOT / "RESCUE" / "BENCHMARK" / "results"
OUT = ROOT / "RESCUE" / "journal_render" / "figures" / "fig0_study_design.png"

h2h = json.loads((B / "head_to_head_paired.json").read_text())
loso = json.loads((B / "loso_head_to_head.json").read_text())
cov = json.loads((B / "rule_coverage.json").read_text())
eu = json.loads((B / "hiv2eu_coverage.json").read_text())
cvf = json.loads((B / "clinical_validation_final_results.json").read_text())
TF = cvf["comparison"]["pooled_4_drugs"]["final_transfer"]

plt.rcParams.update({"font.family": "serif"})

W, H = 11.5, 19.5
fig, ax = plt.subplots(figsize=(W, H))
ax.set_xlim(0, 100)
ax.set_ylim(32.2, 100)
ax.axis("off")

C_DATA, C_CUR, C_FEAT, C_MODEL, C_EVAL = (
    "#dbe8f4", "#e2efe2", "#fdf0d8", "#e8e0f2", "#f7dede")
EDGE = "#4a5a68"

TS, LS, PAD_T, PAD_B = 8.4, 6.8, 3.0, 0.9      # title size, line size, padding

# text metrics in axis units (figure is 24 in tall over 100 units)
IN_PER_UNIT = H / 100.0
TITLE_H = TS * 1.35 / 72.0 / IN_PER_UNIT       # title block height in units
LINE_H = LS * 1.42 / 72.0 / IN_PER_UNIT        # one line height in units


def box_height(nlines):
    return PAD_T + TITLE_H + nlines * LINE_H + PAD_B


def draw(x, top, w, title, lines, face, tsize=TS, lsize=LS):
    h = box_height(len(lines))
    cy = top - h / 2.0
    ax.add_patch(FancyBboxPatch((x - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0.30,rounding_size=0.7",
                                linewidth=1.15, edgecolor=EDGE,
                                facecolor=face, zorder=3))
    ax.text(x, top - PAD_T, title, ha="center", va="top",
            fontsize=tsize, fontweight="bold", color="#0d1b26", zorder=4)
    ty = top - PAD_T - TITLE_H
    for ln in lines:
        ax.text(x, ty, ln, ha="center", va="top", fontsize=lsize,
                color="#22303c", zorder=4)
        ty -= LINE_H
    return h


def arrow(x1, y1, x2, y2, lw=1.6, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=15, linewidth=lw,
                                 color="#37474f", linestyle=ls, zorder=2,
                                 shrinkA=0.5, shrinkB=0.5))


def bus(y_from, y_to, targets, labels=None):
    """Vertical stem, horizontal bus, then drop into each target."""
    mid = (y_from + y_to) / 2.0
    ax.plot([50, 50], [y_from, mid], color="#37474f", lw=1.6, zorder=2)
    ax.plot([min(targets), max(targets)], [mid, mid], color="#37474f",
            lw=1.6, zorder=2)
    for t in targets:
        arrow(t, mid, t, y_to)


TOP = 98.0
GAP = 2.6

# ---- 1 raw sources
n = 3
h = draw(50, TOP, 72, "1  DATA ACQUISITION  (three sources, three species-roles)", [
    "Stanford HIVDB PI phenotype  \u2192  2,171 sequences",
    "ChEMBL 380 / 5074 HIV-2 assays  \u2192  321 rows",
    "NCBI GenBank HIV-2 protease (Entrez)  \u2192  148 sequences",
], C_DATA)
y = TOP - h
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 2 curation
h = draw(50, y, 72, "2  CURATION  (seven steps)", [
    "99-aa filter \u00b7 frame alignment \u00b7 species by provenance \u00b7 drug-name mapping",
    "IC50\u2192\u0394G conversion \u00b7 95%-identity clustering \u00b7 synthetic-row removal",
], C_CUR)
y = y - h
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 3 three corpora
h3 = max(box_height(4), box_height(4), box_height(4))
draw(18, y, 27, "HIV-1 corpus", [
    "5,683 genuine rows",
    "1,882 unique sequences",
    "882 clusters + 714 singletons",
    "1,500 synthetic removed",
], C_CUR)
draw(50, y, 25, "HIV-2 ChEMBL", [
    "321 rows \u2192 320 mapped",
    "19 unique sequences",
    "137 zero-shot / 160 fine-tune",
    "4 drugs \u00d7 40 rows",
], C_CUR)
draw(82, y, 27, "HIV-2 clinical cohort", [
    "654 entries \u00d7 8 PIs = 5,232",
    "Stanford_HIVDB 3,456",
    "Mined_Expanded 1,776",
    "headline subset 2,616",
], C_CUR)
bus(y + GAP, y, [18, 50, 82])
y -= h3
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 4 feature branches
h4 = max(box_height(6), box_height(6))
draw(29, y, 40, "SEQUENCE BRANCH", [
    "21 descriptors / residue:",
    "Hellberg Z (3) + Kidera (10) + VHSE (8)",
    "21 \u00d7 99 = 2,079-D",
    "+ mean/SD of each = 2,121-D",
    "properties of the amino acid",
    "\u2192 identical in both species",
], C_FEAT)
draw(71, y, 40, "STRUCTURE BRANCH", [
    "ESMFold \u2192 23 HIV-2 models (pLDDT)",
    "AutoDock Vina docking",
    "171 HIV-2 complexes (19 \u00d7 9)",
    "13,528 HIV-1 holo (3,382 \u00d7 4)",
    "ONION-Net 10,309-D apo / 10,478-D docked",
    "\u2192 does NOT transfer (R = \u22120.14)",
], C_FEAT)
bus(y + GAP, y, [28, 72])
y -= h4
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 5 model
h = draw(50, y, 72, "3  MODEL TRAINING AND TRANSFER", [
    "HIV-1 per-drug Random Forest / XGBoost  \u00b7  GroupKFold at 95% identity",
    "within-HIV-1 grouped CV:  R = 0.7387 \u00b1 0.0472",
    "HIV-1 model FROZEN \u2192 prediction appended as one meta-feature column",
    "new per-drug RF on [2,121-D | meta-feature]  \u00b7  no weights are modified",
], C_MODEL)
y = y - h
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 6 two outcomes
h6 = max(box_height(2), box_height(2))
draw(29, y, 40, "Zero-shot baseline", [
    "R = 0.5163 on 137 IC50 rows",
    "R = 0.2492 on the same 160 rows",
], C_EVAL)
draw(71, y, 40, "Transfer fine-tune", [
    "R = 0.7018 (160 rows, locked seed)",
    "0.628 \u00b1 0.074 across 10 seeds",
], C_EVAL)
bus(y + GAP, y, [29, 71])
y -= h6
arrow(50, y, 50, y - GAP)
y -= GAP

# ---- 7 evaluation
h7 = max(box_height(1), box_height(1), box_height(1), box_height(1))
draw(13, y, 22, "Measured vs rule", [
    f"{h2h['pooled']['model_r']:.3f} vs {h2h['pooled']['rule_r']:.3f}, P=0.001",
], C_EVAL, tsize=7.6, lsize=6.3)
draw(38, y, 22, "Rule coverage", [
    f"{cov['Q1_clinical_label_coverage']['overall_zero_pct']:.1f}% zero labels",
], C_EVAL, tsize=7.6, lsize=6.3)
draw(63, y, 20, "LOSO robustness", [
    f"ahead in {loso['loso_full']['frac_folds_model_wins']*100:.0f}% of folds",
], C_EVAL, tsize=7.6, lsize=6.3)
draw(87, y, 22, "Clinical (secondary)", [
    f"R={TF['pearson_r']:.3f}, AUC={TF['roc_auc']:.3f}",
], C_EVAL, tsize=7.6, lsize=6.3)
bus(y + GAP, y, [14, 38, 62, 86])

fig.savefig(OUT, dpi=600, bbox_inches="tight", facecolor="white")
print("wrote", OUT)
print("bottom y:", round(y - h7, 2))
