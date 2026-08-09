#!/usr/bin/env python3
"""Generate the two-panel kinship summary figure.
Panel A: TKGWV2 / KING-robust / READv2 reported coefficients and method-specific thresholds.
Panel B: KING kinship-vs-IBS0 scatter placing ORT15xORT16 at parent-offspring.
"""
import argparse
import csv
from pathlib import Path
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser()
parser.add_argument(
    "--tkgwv2", type=Path,
    default=ROOT / "data/summary/kinship/tkgwv2_pair_results.tsv",
)
parser.add_argument(
    "--readv2", type=Path,
    default=ROOT / "data/summary/kinship/readv2_ort15_ort16.tsv",
)
parser.add_argument(
    "--king-ibs0", type=Path,
    default=ROOT / "data/summary/kinship/king_ibs0_summary.tsv",
)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

def table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


tkgwv2 = next(
    row for row in table(args.tkgwv2)
    if row["sample1"] == "ORT15_fu_bamrefine5" and row["sample2"] == "ORT16_fu_bamrefine5"
)
readv2 = next(row for row in table(args.readv2) if row["analysis"] == "primary_chr1_22_X_Y")
king = next(row for row in table(args.king_ibs0) if row["analysis"] == "genomewide")

highlight_color = "#c1121f"
methods=["TKGWV2","KING-robust","READv2"]
coeffs=[
    float(tkgwv2["HRC"]),
    float(king["KING_robust"]),
    round(float(readv2["KinshipCoefficient"]), 3),
]
snps=[
    f"{int(tkgwv2['used_snps']) / 1e6:.1f}M SNPs",
    f"{int(king['sites']) / 1e6:.1f}M SNPs",
    f"{round(int(readv2['OverlapNSNPs']) / 1000):d}k sites",
]
king_kin=float(king["KING_robust"])
ibs0_obs=float(king["observed_IBS0_per_site"])
ibs0_fs=float(king["full_sibling_expected"])
ibs0_unrel=float(king["unrelated_expected"])
ref_pts={"MZ twin / duplicate":(0.50,0.0),"Parent-offspring":(0.25,0.0),
         "Full siblings":(0.25,ibs0_fs),"2nd degree":(0.125,ibs0_fs*2),"Unrelated":(0.0,ibs0_unrel)}
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.6,3.7),gridspec_kw={"width_ratios":[1.0,1.15]})
axA.axhline(0.177,color="#238b45",lw=1.0,ls="--")
axA.axhline(0.1875,color="#2b6cb0",lw=1.0,ls=":")
axA.text(-0.52,0.171,"KING boundary 0.177",ha="left",va="top",fontsize=5.8,color="#238b45")
axA.text(-0.52,0.193,"TKGWV2 boundary 0.188",ha="left",va="bottom",fontsize=5.8,color="#2b6cb0")
x=np.arange(3); axA.scatter(x,coeffs,s=70,color=highlight_color,zorder=5,edgecolor="black",linewidth=0.5)
for xi,c,s in zip(x,coeffs,snps):
    axA.annotate(f"{c:.3f}",(xi,c),textcoords="offset points",xytext=(0,8),ha="center",fontsize=6.5,fontweight="bold")
    axA.annotate(s,(xi,c),textcoords="offset points",xytext=(0,-13),ha="center",fontsize=5.5,color="#555555")
axA.set_xticks(x); axA.set_xticklabels(methods, fontsize=7.5); axA.set_ylim(0,0.40); axA.set_xlim(-0.6,2.6)
axA.set_ylabel("Reported relatedness coefficient",labelpad=10); axA.set_title("Three estimates fall in the first-degree range",loc="left",fontsize=8)
axA.annotate("a",xy=(0,1),xycoords="axes fraction",xytext=(-42,10),textcoords="offset points",fontsize=11,fontweight="bold",va="bottom",annotation_clip=False)
for name,(kx,ky) in ref_pts.items():
    axB.scatter(kx,ky,s=45,facecolor="white",edgecolor="#888888",linewidth=0.9,zorder=3)
axB.scatter(king_kin,ibs0_obs,s=95,color=highlight_color,edgecolor="black",linewidth=0.6,zorder=6)
axB.annotate("ORT15 × ORT16\n(0.238, 0.0034)",(king_kin,ibs0_obs),textcoords="offset points",xytext=(-6,10),ha="right",fontsize=6.5,color=highlight_color,fontweight="bold")
lab_off={"MZ twin / duplicate":(0,9,"center","bottom"),"Parent-offspring":(6,6,"left","center"),"Full siblings":(8,0,"left","center"),"2nd degree":(8,0,"left","center"),"Unrelated":(6,0,"left","center")}
for name,(kx,ky) in ref_pts.items():
    dx,dy,ha,va=lab_off[name]; axB.annotate(name,(kx,ky),textcoords="offset points",xytext=(dx,dy),ha=ha,va=va,fontsize=6.8,color="#555555")
axB.set_xlim(-0.03,0.62); axB.set_ylim(-0.004,0.082)
axB.set_xlabel("KING kinship coefficient"); axB.set_ylabel("IBS0 rate (per site)",labelpad=8)
axB.set_title("IBS0 favours parent-offspring",loc="left",fontsize=8)
axB.annotate("b",xy=(0,1),xycoords="axes fraction",xytext=(-42,10),textcoords="offset points",fontsize=11,fontweight="bold",va="bottom",annotation_clip=False)
fig.tight_layout()
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, bbox_inches="tight")
print(f"wrote {args.output}")
