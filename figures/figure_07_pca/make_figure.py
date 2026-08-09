#!/usr/bin/env python3
"""Generate the three-panel PCA figure from compact aggregate inputs.

Datasets: 24 unique treatment datasets total (12 per individual: four full-UDG +
eight non-UDG), each in two representations (24 imputed + 24 pseudo-haploid) = 48
projected points on one shared reference basis.

Reference backgrounds are population-level aggregate density masks and fixed-bin
counts rather than individual AADR coordinates. Orientation signs are cosmetic;
distances and explained-variance percentages are unchanged.
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300, "font.size": 8,
    "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 7, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.8, "font.family": "sans-serif",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
QUERY = os.path.join(HERE, "inputs", "query_coordinates.tsv")
EVAL = os.path.join(ROOT, "data", "summary", "pca", "pca_aadr_matched48.eval")
AGGREGATES = os.path.join(HERE, "inputs", "reference_aggregates.npz")
OUTPUT = os.environ.get("ORT_FIGURE_OUTPUT", os.path.join(HERE, "figure_07_pca.pdf"))

SX, SY = -1.0, -1.0  # orientation flip (cosmetic)

def load_queries(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append((row["iid"], [float(row["pc1"]), float(row["pc2"])],
                         row["representation"]))
    return rows

rows = load_queries(QUERY)
aggregate = np.load(AGGREGATES, allow_pickle=False)
evals = [float(x) for x in open(EVAL).read().split()]
trace = sum(v for v in evals if v > 0)
var_pct = [100.0 * evals[i] / trace for i in range(len(evals))]
PC1LAB = f"PC1 ({var_pct[0]:.2f}%)"
PC2LAB = f"PC2 ({var_pct[1]:.2f}%)"

imp  = [(iid, pcs) for iid, pcs, representation in rows if representation == "IMP"]
phap = [(iid, pcs) for iid, pcs, representation in rows if representation == "PHAP"]
print(f"query rows: imp={len(imp)} phap={len(phap)} total={len(rows)}")
assert len(imp) == 24 and len(phap) == 24

def X(pcs): return SX * pcs[0]
def Y(pcs): return SY * pcs[1]

# ---- Shared continent palette (saturated Set1) ----
CONT_ORDER = ["Africa","WestEurasia","SouthAsia","CentralAsiaSiberia","EastAsia","Oceania","America"]
CONT_COL = {"Africa":"#a65628","WestEurasia":"#e41a1c","SouthAsia":"#4daf4a",
            "CentralAsiaSiberia":"#377eb8","EastAsia":"#984ea3","Oceania":"#f781bf",
            "America":"#ff7f00"}
TREAT_ORDER = ["raw","trim5","trim10","rescale5","rescale10","rescaled","bamrefine5","bamrefine10"]
TREAT_COLORS = {"raw":"#000000","trim5":"#1f77b4","trim10":"#3690c0","rescale5":"#2ca02c",
                "rescale10":"#74c476","rescaled":"#006d2c","bamrefine5":"#d62728","bamrefine10":"#ff7f0e"}
def parse_query(iid):
    base = iid.split("__", 1)[1]; parts = base.split("_")
    return parts[0], parts[1], "_".join(parts[2:])
def marker_for(indiv):
    return "o" if indiv == "ORT15" else "s"

conts = CONT_ORDER
YAK = tuple(float(value) for value in aggregate["yakut_mean"])
YAK_N = int(aggregate["yakut_n"][0])
zoom_window = aggregate["zoom_window"]
ZX = (float(zoom_window[0]), float(zoom_window[1]))
ZY = (float(zoom_window[2]), float(zoom_window[3]))

def draw_aggregate(ax, panel, continent, alpha, dot_size):
    """Draw only disclosure-controlled aggregate reference information."""
    prefix = f"{panel}_{continent}"
    xs = aggregate[f"{prefix}_x"]
    ys = aggregate[f"{prefix}_y"]
    mask = aggregate[f"{prefix}_mask"].astype(bool)
    if mask.any():
        xx, yy = np.meshgrid(xs, ys)
        ax.contourf(xx, yy, mask.astype(float), levels=[0.5, 1.5],
                    colors=[CONT_COL[continent]], alpha=alpha, zorder=1)
    bins = aggregate[f"{prefix}_bins"]
    if len(bins):
        sizes = dot_size + 1.8 * np.sqrt(bins[:, 2])
        ax.scatter(bins[:, 0], bins[:, 1], s=sizes, c=CONT_COL[continent],
                   edgecolors="none", alpha=0.70, rasterized=True, zorder=2)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(6.5, 6.2))
gs = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.0], hspace=0.42, wspace=0.30,
                      left=0.115, right=0.985, top=0.955, bottom=0.175)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])

# ---- Panel A: worldwide reference (aggregate colour plates and dots) ----
for c in conts:
    draw_aggregate(axA, "world", c, alpha=0.12, dot_size=2.8)
axA.add_patch(Rectangle((ZX[0],ZY[0]), ZX[1]-ZX[0], ZY[1]-ZY[0], fill=False,
              edgecolor="black", lw=1.0, ls="--", zorder=6))
# two gold stars = pseudo-haploid full-UDG raw anchor per individual
cy = float(np.mean([Y(pcs) for _, pcs in imp + phap]))
for iid, pcs in phap:
    ind, lib, tr = parse_query(iid)
    if lib == "fu" and tr == "raw":
        axA.scatter(X(pcs), Y(pcs), s=95, marker="*", facecolors="#ffd11a",
                    edgecolors="black", linewidths=0.5, zorder=7)
axA.annotate("ORT15/ORT16\n(projected)", xy=(ZX[0], cy),
             xytext=(ZX[0]-0.9*(ZX[1]-ZX[0])-0.006, cy+0.006),
             fontsize=6.6, ha="right", va="center",
             arrowprops=dict(arrowstyle="->", lw=0.7, color="black"), zorder=8)
axA.set_xlabel(PC1LAB); axA.set_ylabel(PC2LAB)
axA.set_title("AADR v62 worldwide reference", loc="left")
cont_h = [Line2D([0],[0], marker='o', ls='none', mfc=CONT_COL[c], mec='none', ms=4.5, label=c)
          for c in conts]
cont_h.append(Line2D([0],[0], marker='*', ls='none', mfc='#ffd11a', mec='black', mew=0.4,
              ms=8, label='ORT15/ORT16'))
axA.legend(handles=cont_h, fontsize=5.8, loc="upper left", frameon=True, framealpha=0.8,
           edgecolor="none", handletextpad=0.2, labelspacing=0.22, borderpad=0.3)
axA.margins(0.05)

# ---- Panels B & C: zooms with aggregate colour plates, reference dots, and ORT ----
def draw_zoom(ax, data, title):
    ext = (ZX[0], ZX[1], ZY[0], ZY[1])
    for c in conts:
        draw_aggregate(ax, "zoom", c, alpha=0.16, dot_size=5.0)
    ax.scatter(*YAK, marker="X", s=64, c="#111111", edgecolors="white", linewidths=0.6, zorder=3)
    for iid, pcs in data:
        indiv, lib, treat = parse_query(iid)
        col = TREAT_COLORS.get(treat, "#555555"); mk = marker_for(indiv)
        if lib == "fu":
            ax.scatter(X(pcs), Y(pcs), s=48, marker=mk, facecolors=col,
                       edgecolors="black", linewidths=0.35, alpha=0.95, zorder=4)
        else:
            ax.scatter(X(pcs), Y(pcs), s=48, marker=mk, facecolors="none",
                       edgecolors=col, linewidths=1.4, zorder=4)
    ax.set_xlim(*ZX); ax.set_ylim(*ZY)
    ax.set_xlabel(PC1LAB); ax.set_ylabel(PC2LAB); ax.set_title(title, loc="left")
    # plain decimal ticks (no x10^-2 offset text)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(4))
    ax.yaxis.set_major_locator(mticker.MaxNLocator(4))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
    ax.tick_params(labelsize=6.6)
draw_zoom(axB, imp,  "Imputed (24 datasets)")
draw_zoom(axC, phap, "Pseudo-haploid (24 datasets)")

# ---- shared bottom legend (treatments + ORT/library encoding + Yakut) ----
# "rescaled" is mapDamage --rescale with the default 12-bp window.
TREAT_LABEL = {
    "raw": "Uncorrected",
    "trim5": "Trim-5",
    "trim10": "Trim-10",
    "rescale5": "Rescale-5",
    "rescale10": "Rescale-10",
    "rescaled": "Rescale-12",
    "bamrefine5": "bamRefine-5",
    "bamrefine10": "bamRefine-10",
}
leg_t = [Line2D([0],[0], marker='s', ls='none', mfc=TREAT_COLORS[t], mec='black', mew=0.3,
         ms=6, label=TREAT_LABEL.get(t, t)) for t in TREAT_ORDER]
leg_e = [Line2D([0],[0], marker='o', ls='none', mfc='none', mec='black', ms=6.5, label='ORT15'),
         Line2D([0],[0], marker='s', ls='none', mfc='none', mec='black', ms=6.5, label='ORT16'),
         Line2D([0],[0], marker='o', ls='none', mfc='dimgray', mec='black', ms=6.5, label='full-UDG (filled)'),
         Line2D([0],[0], marker='o', ls='none', mfc='none', mec='black', ms=6.5, label='non-UDG (hollow)'),
         Line2D([0],[0], marker='X', ls='none', mfc='#111111', mec='white', mew=0.5,
                ms=7, label=f'Yakut mean (n={YAK_N})')]
fig.legend(handles=leg_t+leg_e, fontsize=6.6, loc="lower center", ncol=5,
           bbox_to_anchor=(0.5, 0.0), frameon=False, handletextpad=0.3,
           columnspacing=1.0, labelspacing=0.35)

for ax, L in ((axA,'a'), (axB,'b'), (axC,'c')):
    ax.annotate(L, xy=(0,1), xycoords="axes fraction", xytext=(-32,4),
                textcoords="offset points", fontsize=11, fontweight="bold", va="bottom")

os.makedirs(os.path.dirname(os.path.abspath(OUTPUT)), exist_ok=True)
fig.savefig(OUTPUT)
print(f"wrote {OUTPUT}")
print(f"PC1={var_pct[0]:.4f}% PC2={var_pct[1]:.4f}%  flip=({SX},{SY})  figsize=(6.5,6.2)")
print(f"Yakut mean={YAK}")
