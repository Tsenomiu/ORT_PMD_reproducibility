#!/usr/bin/env python3
"""Plot imputation concordance by MAF bin with an all-site NRD inset.

Dosage r-squared uses column c5 of the per-bin GLIMPSE summary; c4 is the
best-guess statistic. NRD bars begin at zero. The plotted NRD values are checked
against an independently harvested GCsV summary for all 16 sample-treatment rows.
"""
import argparse, json, csv, os
import matplotlib; matplotlib.use("Agg")
# Embed TrueType (Type 42) fonts, not Matplotlib Type 3 — glyphs identical, format only.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes as _inset_axes
HERE=os.path.dirname(os.path.abspath(__file__))

rs=json.load(open(os.path.join(HERE,"rsquare_all16.json")))
TAB={}
# Independent concordance GCsV harvest (NRD = field 0).
GCSV=json.load(open(os.path.join(HERE,"impute_gcsv.json")))  # {ind:{treatment:(NRD,best-guess,dosage)}}

# (display, key, colour, ls, marker) — original Fig 5 styling
STRATEGIES=[
    ("Raw",          "raw",         "#7f7f7f", ":",  "x"),
    ("Trim-5",       "trim5",       "#ff7f0e", "-",  "o"),
    ("Trim-10",      "trim10",      "#d62728", "--", "s"),
    ("Rescale-5",    "rescale5",    "#1f77b4", "-",  "o"),
    ("Rescale-10",   "rescale10",   "#17becf", "--", "s"),
    ("Rescale-12",   "rescaled",    "#9467bd", ":",  "^"),
    ("bamRefine-5",  "bamrefine5",  "#2ca02c", "-",  "D"),
    ("bamRefine-10", "bamrefine10", "#bcbd22", "--", "P"),
]
BIN_ORDER=list(range(1,9))

def curve(sample,key):
    b=[x for x in rs[f"{sample}_{key}"] if x["bin"]>=1]
    return [x["meanMAF"] for x in b],[x["c5"] for x in b]   # c5 = DOSAGE

def make_inset(ax_main, nrd_data):
    ax_in=_inset_axes(ax_main,width="46%",height="40%",loc="lower right",
                      bbox_to_anchor=(0.0,0.12,1.0,1.0),bbox_transform=ax_main.transAxes,borderpad=0)
    names=[n for n,_ in nrd_data]; nrds=[v for _,v in nrd_data]
    colours=[s[2] for s in STRATEGIES]
    x=np.arange(len(names))
    bars=ax_in.bar(x,nrds,color=colours,edgecolor="black",lw=0.4,width=0.7)
    ax_in.set_xticks(x)
    ax_in.set_xticklabels([n.replace("bamRefine","bR").replace("Rescale","R").replace("Trim","T") for n in names],
                          rotation=90,ha="center",fontsize=5)
    ax_in.set_ylabel("NRD (%)",fontsize=10); ax_in.set_title("NRD (all sites)",fontsize=10,pad=2)
    ax_in.tick_params(axis="both",labelsize=9); ax_in.tick_params(axis="x",pad=2)
    ax_in.grid(axis="y",ls=":",alpha=0.4)
    for b,v in zip(bars,nrds):
        ax_in.text(b.get_x()+b.get_width()/2, b.get_height()+0.15, f"{v:.1f}",
                   ha="center",va="bottom",fontsize=5)
    ax_in.set_ylim(0, max(nrds)*1.15)   # BARS FROM ZERO (no truncated axis)

def plot_panel(ax, sample, title, show_legend):
    first_mafs,_=curve(sample,STRATEGIES[0][1])
    nrd_data=[]
    for name,key,colour,ls,marker in STRATEGIES:
        mafs,r2s=curve(sample,key)
        ax.plot(mafs,r2s,color=colour,ls=ls,marker=marker,ms=5,lw=1.5,label=name)
        nrd=round(float(TAB[(sample,key)]["NRD_pct"]),2)
        nrd_data.append((name,nrd))
    ax.set_xscale("log"); ax.set_xlim(0.002,0.6); ax.set_ylim(0.0,1.02)
    ax.set_xlabel("Minor allele frequency (MAF)",fontsize=14)
    ax.set_ylabel("Dosage $r^2$",fontsize=14); ax.set_title(title,fontsize=12,pad=18)
    ax.grid(ls=":",alpha=0.35)
    ax2=ax.twiny(); ax2.set_xscale("log"); ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(first_mafs); ax2.set_xticklabels([f"bin{b}" for b in BIN_ORDER],fontsize=8,rotation=90,ha="center")
    ax2.tick_params(axis="x",length=3,pad=2)
    if show_legend: ax.legend(fontsize=11,loc="upper left",framealpha=0.9)
    make_inset(ax,nrd_data)
    return nrd_data

def main():
    global TAB
    parser=argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    args=parser.parse_args()
    recs=list(csv.DictReader(open(args.summary)))
    TAB={(r["sample"],r["treatment"]):r for r in recs}
    output=os.path.abspath(args.output)
    os.makedirs(os.path.dirname(output),exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(13,5.5))
    fig.subplots_adjust(wspace=0.35,top=0.84,bottom=0.12,left=0.07,right=0.97)
    nd15=plot_panel(axes[0],"ORT15","ORT15 non-UDG (0.27\u00d7) vs matched full-UDG comparator (0.90\u00d7)",True)
    nd16=plot_panel(axes[1],"ORT16","ORT16 non-UDG (0.35\u00d7) vs matched full-UDG comparator (1.46\u00d7)",False)

    # ---- assertions: every plotted NRD (from TAB) must match the INDEPENDENT GCsV
    #      harvest to 2 dp, so a bad table value fails generation (req 4). ----
    name2key=dict((s[0],s[1]) for s in STRATEGIES)
    n=0
    for sample,nd in [("ORT15",nd15),("ORT16",nd16)]:
        for name,nrd in nd:                       # nrd is the plotted value (read from TAB)
            key=name2key[name]
            src=round(float(GCSV[sample][key][0]),2)   # independent source of truth
            assert nrd==src, f"NRD mismatch {sample}/{name}: plotted/table {nrd} vs GCsV {src}"
            n+=1
    assert n==16, n

    fig.savefig(output,dpi=300,bbox_inches="tight")

    plt.close(fig)

    # ---- per-bin plotting-data TSV ----
    rows=[]
    for sample in ["ORT15","ORT16"]:
        for _,key,_,_,_ in STRATEGIES:
            for x in rs[f"{sample}_{key}"]:
                if x["bin"]>=1:
                    rows.append((sample,key,x["bin"],round(x["meanMAF"],6),round(x["c5"],6),round(x["c4"],6)))
    with open(os.path.join(os.path.dirname(output),"figure_05_plotdata_perbin.tsv"),"w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t")
        w.writerow(["sample","treatment","maf_bin","mean_MAF","dosage_r2","bestguess_r2"]); w.writerows(rows)

    print("=== Figure 5 fixed-comparator verification ===")
    print(f"NRD plotted == table CSV (2dp) for all {n}: PASS")
    for sample in ["ORT15","ORT16"]:
        sub=[(nm,float(TAB[(sample,k)]["dosage_r2_maf01"]),float(TAB[(sample,k)]["NRD_pct"])) for nm,k,_,_,_ in STRATEGIES]
        bestds=max(sub,key=lambda z:z[1]); lown=min(sub,key=lambda z:z[2])
        print(f"{sample}: best dosage r2 = {bestds[0]} {bestds[1]:.4f} | lowest NRD = {lown[0]} {lown[2]:.2f}")
    print(f"per-bin plotdata rows: {len(rows)} (expect 128); PDF and TSV written.")

if __name__=="__main__": main()
