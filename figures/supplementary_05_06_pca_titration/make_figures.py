#!/usr/bin/env python3
"""Generate the ORT15 and ORT16 coverage-titration PCA figures.

Deterministic regeneration from disclosure-controlled derived outputs of a
SEPARATE smartpca run from the matched-48 main PCA. Both library types (full-UDG
and non-UDG) of each individual were
downsampled to target coverages 0.01-0.2x (10 replicates per reachable target) plus
the full-coverage anchor, each PMD treatment pseudo-haploid-called (pileupCaller
--randomHaploid, 1240k) and projected onto the AADR v62 HGDP+SGDP reference PCA.

Each individual is split into two stacked coverage sections: 0.01–0.06x and
0.08x–full. Each section contains full-UDG and non-UDG rows across four coverage
columns.

ORT coordinates, orientation (ROT=-1 negate both PCs), continent palette,
treatment colours, aggregate continent plates and the tight ORT zoom window are
unchanged. Reference-individual identifiers and coordinates are not distributed;
the inputs contain density masks and fixed-bin counts with a minimum cell count of
three.

Inputs: inputs/query_coordinates.tsv, inputs/reference_aggregates.npz
Outputs: supplementary_05_pca_titration_ort15.pdf and
         supplementary_06_pca_titration_ort16.pdf
The plotting code is compatible with the locally available Matplotlib stack.
The PDF wraps a 600-dpi PNG with deterministic pdfTeX settings because the
Matplotlib PDF backend clips one repeated-marker panel.
"""
import csv, os, re, shutil, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)) or "."
INPUTS = os.path.join(HERE, "inputs")
QUERY = os.path.join(INPUTS, "query_coordinates.tsv")
AGGREGATES = os.path.join(INPUTS, "reference_aggregates.npz")
OUTPUT_DIR = os.environ.get("ORT_FIGURE_OUTPUT_DIR", HERE)
PDFTEX = shutil.which("pdftex")

# Two coverage sections; one supplementary figure per individual.
SECTIONS = [("0.01--0.06$\\times$", ["0.01","0.02","0.04","0.06"]),
            ("0.08$\\times$--full", ["0.08","0.1","0.2","full"])]
METH_ORDER = ["raw","trim5","trim10","rescale5","rescale10","rescaled","bamrefine5","bamrefine10"]
METH_COL = {"raw":"#000000","trim5":"#1f77b4","trim10":"#3690c0","rescale5":"#2ca02c",
            "rescale10":"#74c476","rescaled":"#006d2c","bamrefine5":"#d62728","bamrefine10":"#ff7f0e"}
LIBNAME = {"fu":"full-UDG","nu":"non-UDG"}
CONT_ORDER = ["Africa","WestEurasia","SouthAsia","CentralAsiaSiberia","EastAsia","Oceania","America"]
CONT_COL = {"Africa":"#a65628","WestEurasia":"#e41a1c","SouthAsia":"#4daf4a",
            "CentralAsiaSiberia":"#377eb8","EastAsia":"#984ea3","Oceania":"#f781bf","America":"#ff7f00"}
ROT = -1.0   # negate both PCs (conventional orientation) — coordinates unchanged in magnitude

matplotlib.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 300, "font.size": 8, "font.family": "sans-serif",
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.6,
    "svg.fonttype": "path", "svg.hashsalt": "ORT-PCA-titration",
})

def parse(iid):
    m = re.match(r"(ORT1[56])_(fu|nu)_(\w+)__c([0-9.]+|full)__r(\d+)$", iid)
    if not m: return None
    return dict(ind=m.group(1), lib=m.group(2), meth=m.group(3), cov=m.group(4), rep=int(m.group(5)))

def load():
    ort = []
    with open(QUERY, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = parse(row["iid"])
            if sample:
                ort.append((sample, ROT*float(row["pc1"]), ROT*float(row["pc2"])))
    aggregate = np.load(AGGREGATES, allow_pickle=False)
    yak_mean = tuple(float(value) for value in aggregate["yakut_mean"])
    n_yak = int(aggregate["yakut_n"][0])
    return ort, aggregate, yak_mean, n_yak

def draw_cont(ax, aggregate, ind, dot_size=5, alpha=0.5, hull_alpha=0.10, visible=None):
    """Draw disclosure-controlled aggregate reference density and bin counts."""
    for c in CONT_ORDER:
        prefix = f"{ind}_{c}"
        xs = aggregate[f"{prefix}_x"]
        ys = aggregate[f"{prefix}_y"]
        mask = aggregate[f"{prefix}_mask"].astype(bool)
        if hull_alpha > 0 and mask.any():
            xx, yy = np.meshgrid(xs, ys)
            ax.contourf(xx, yy, mask.astype(float), levels=[0.5, 1.5],
                        colors=[CONT_COL[c]], alpha=hull_alpha, zorder=1)
        bins = aggregate[f"{prefix}_bins"]
        if len(bins):
            sizes = dot_size + 1.8 * np.sqrt(bins[:, 2])
            ax.scatter(bins[:, 0], bins[:, 1], s=sizes, c=CONT_COL[c],
                       lw=0, alpha=alpha, zorder=2, rasterized=True)
        if visible is not None and bool(aggregate[f"{prefix}_visible"][0]):
            visible.add(c)

def ortwin(ort, pf=0.18, pa=1e-3):
    xs=[x for _,x,_ in ort]; ys=[y for _,_,y in ort]
    px=(max(xs)-min(xs))*pf+pa; py=(max(ys)-min(ys))*pf+pa
    return (min(xs)-px, max(xs)+px, min(ys)-py, max(ys)+py)

def build(ind, ort, aggregate, yak_mean, n_yak, out_base):
    win = tuple(float(value) for value in aggregate[f"{ind}_window"])
    calculated_win = ortwin(ort)
    assert np.allclose(win, calculated_win, atol=1e-12), (win, calculated_win)
    # Explicit four-sided polygon keeps the ORT16 symbol square while avoiding
    # a Matplotlib SVG/PDF marker-id collision triggered by marker="s".
    marker = "o" if ind == "ORT15" else (4, 0, 45)
    # two stacked sections; each = 2 rows (fu/nu) x 4 cols (coverage)
    # authored at final display size (6.5in wide = \linewidth) so LaTeX applies NO
    # downscaling and in-figure point sizes render 1:1 on the page
    fig = plt.figure(figsize=(6.5, 6.15))
    outer = fig.add_gridspec(2, 1, hspace=0.40, left=0.14, right=0.985, top=0.912, bottom=0.255)
    present_meth = set()
    visible_cont = set()
    block_axes = {}
    for si,(sect_title, covs) in enumerate(SECTIONS):
        gs = outer[si].subgridspec(2, len(covs), hspace=0.10, wspace=0.10)
        block_axes[si] = []
        for ri, lib in enumerate(("fu","nu")):
            for ci, cov in enumerate(covs):
                ax = fig.add_subplot(gs[ri, ci])
                # Rasterise only the descriptive reference layer (zorder < 3).
                # Axes, labels and all ORT treatment markers remain vector.
                ax.set_rasterization_zorder(3)
                block_axes[si].append(ax)
                draw_cont(ax, aggregate, ind, dot_size=3.2, alpha=0.5,
                          hull_alpha=0.10, visible=visible_cont)
                # Yakut mean marker (nearest well-sampled reference pop), as main Figure 9 b/c
                if yak_mean is not None:
                    ax.scatter([yak_mean[0]],[yak_mean[1]], marker="X", s=26, c="black",
                               lw=0.4, edgecolor="white", zorder=5)
                for meth in METH_ORDER:
                    pts=[(x,y) for s,x,y in ort
                         if s["ind"]==ind and s["lib"]==lib and s["meth"]==meth and s["cov"]==cov]
                    if pts:
                        present_meth.add(meth)
                        col = METH_COL.get(meth,"#888")
                        # fill encodes library: filled=full-UDG ("full"), hollow=non-UDG
                        if lib == "fu":
                            ax.scatter([x for x,_ in pts],[y for _,y in pts], s=16, marker=marker,
                                       c=col, lw=0.25, edgecolor="white", alpha=0.9, zorder=6)
                        else:
                            ax.scatter([x for x,_ in pts],[y for _,y in pts], s=16, marker=marker,
                                       facecolors="none", edgecolors=col, lw=0.9, alpha=0.95, zorder=6)
                ax.set_xlim(win[0],win[1]); ax.set_ylim(win[2],win[3])
                ax.tick_params(labelsize=6.5, length=2, pad=1)
                ax.xaxis.set_major_locator(mticker.MaxNLocator(3))
                ax.yaxis.set_major_locator(mticker.MaxNLocator(3))
                # Three decimals are necessary here: at this zoom, two decimals
                # rounds distinct ticks (for example 0.015 and 0.020) to the same label.
                ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
                ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                if ri==0:
                    ax.set_title(f"{cov}$\\times$" if cov!="full" else "full", fontsize=9, pad=3)
                    ax.set_xticklabels([])
                if ci==0:
                    ax.set_ylabel(f"{LIBNAME[lib]}\nPC2", fontsize=8)
                else:
                    ax.set_yticklabels([])
                if ri==1:
                    ax.set_xlabel("PC1", fontsize=8)
                ax.grid(ls=":", alpha=0.3, lw=0.4)
    # section coverage-range labels at the far left, rotated, spanning each block
    # (clear of the centred per-column coverage titles)
    fig.canvas.draw()
    for si,(sect_title, covs) in enumerate(SECTIONS):
        boxes = [a.get_position() for a in block_axes[si]]
        y_ctr = 0.5*(min(b.y0 for b in boxes) + max(b.y1 for b in boxes))
        fig.text(0.028, y_ctr, sect_title.replace("$\\times$","\u00d7"), fontsize=9.5,
                 fontweight="bold", ha="center", va="center", rotation=90)
    # shared legends along the bottom — same encoding as main Figure 9 b/c
    # "rescaled" = mapDamage --rescale default 12-bp window; label rescale12 to
    # use the Rescale-12 display label (data key unchanged).
    METH_LABEL = {"rescaled": "rescale12"}
    mh=[Line2D([0],[0],marker=marker,ls="",ms=6,mec="white",mew=0.4,mfc=METH_COL[m],label=METH_LABEL.get(m,m))
        for m in METH_ORDER if m in present_meth]
    # library type = fill only (Yakut is a reference population, so it belongs
    # with the modern-population key below, not here). Fill state is obvious from
    # the marker glyph, so labels are kept short (filled=full-UDG, hollow=non-UDG).
    enc=[Line2D([0],[0],marker=marker,ls="",ms=6.5,mfc="#555",mec="white",mew=0.4,label=f"full-UDG"),
         Line2D([0],[0],marker=marker,ls="",ms=6.5,mfc="none",mec="#555",mew=1.0,label=f"non-UDG")]
    ch=[Line2D([0],[0],marker="o",ls="",ms=6,mfc=CONT_COL[c],mec="none",label=c) for c in CONT_ORDER if c in visible_cont]
    ch.append(Line2D([0],[0],marker="X",ls="",ms=6.5,mfc="black",mec="white",mew=0.4,label=f"Yakut mean (n={n_yak})"))
    # library-type fill encodings are folded into the PMD-method legend (no
    # separate "library type" titled block); balanced 2-row grid so the fill
    # entries sit in the last column without overflow.
    allmh = mh+enc
    ncol_meth = -(-len(allmh)//2)  # ceil to 2 rows
    leg1=fig.legend(handles=allmh, loc="lower center", ncol=ncol_meth, fontsize=6.8,
                    title="PMD correction method (colour)",
                    title_fontsize=7.4, frameon=False,
                    bbox_to_anchor=(0.5,0.10), handletextpad=0.25, columnspacing=0.9)
    fig.add_artist(leg1)
    fig.legend(handles=ch, loc="lower center", ncol=len(ch), fontsize=6.8,
               title="Modern population", title_fontsize=7.4, frameon=False,
               bbox_to_anchor=(0.5,0.006), handletextpad=0.25, columnspacing=0.9)
    fig.suptitle(f"{ind}: coverage-titration PCA (both library types) — same encoding as main PCA panels b/c",
                 fontsize=9.5, y=0.965)
    fig.savefig(out_base+".png", dpi=600,
                metadata={"Software": "make_figures.py"})
    # Editable vector source.
    fig.savefig(out_base+".svg", metadata={
        "Title": f"{ind} coverage-titration PCA",
        "Creator": "make_figures.py",
    })
    if PDFTEX is None:
        raise RuntimeError("pdftex is required to create the deterministic PDF wrapper")
    # The available Matplotlib PDF backend clips one repeated-marker panel.
    # Wrap the verified 600-dpi PNG without resampling. At the native 6.5-inch
    # width this exceeds the journal's 300-dpi requirement for colour figures.
    stem = os.path.basename(out_base)
    wrapper = out_base+".wrap.tex"
    tex = (
        "\\pdfoutput=1\n"
        "\\pdfpagewidth=468bp\n"
        "\\pdfpageheight=442.8bp\n"
        "\\pdfhorigin=0bp\n"
        "\\pdfvorigin=0bp\n"
        "\\pdfinfoomitdate=1\n"
        "\\pdfsuppressptexinfo=15\n"
        "\\pdftrailerid{}\n"
        f"\\setbox0=\\hbox{{\\pdfximage width 468bp height 442.8bp {{{stem}.png}}"
        "\\pdfrefximage\\pdflastximage}\n"
        "\\shipout\\box0\n"
        "\\end\n"
    )
    with open(wrapper, "w", encoding="ascii") as fh:
        fh.write(tex)
    subprocess.run([PDFTEX, "-interaction=nonstopmode", "-halt-on-error",
                    os.path.basename(wrapper)], cwd=os.path.dirname(out_base),
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
    os.replace(out_base+".wrap.pdf", out_base+".pdf")
    for suffix in (".wrap.tex", ".wrap.log"):
        path = out_base+suffix
        if os.path.exists(path):
            os.remove(path)
    plt.close(fig)
    return present_meth

def main():
    ort, aggregate, yak_mean, n_yak = load()
    n_by = {}
    for s,_,_ in ort: n_by[s["ind"]] = n_by.get(s["ind"],0)+1
    print("ORT projected rows:", len(ort), "by individual:", n_by,
          "reference background: disclosure-controlled aggregate",
          "Yakut mean:", yak_mean, "n_yak:", n_yak)
    m15 = build("ORT15", [o for o in ort if o[0]["ind"]=="ORT15"], aggregate, yak_mean, n_yak,
                os.path.join(OUTPUT_DIR,"supplementary_05_pca_titration_ort15"))
    m16 = build("ORT16", [o for o in ort if o[0]["ind"]=="ORT16"], aggregate, yak_mean, n_yak,
                os.path.join(OUTPUT_DIR,"supplementary_06_pca_titration_ort16"))
    print("methods present ORT15:", sorted(m15))
    print("methods present ORT16:", sorted(m16))

if __name__ == "__main__":
    main()
