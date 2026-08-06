#!/usr/bin/env python3
"""Generate the two-panel kinship summary figure.
Panel A: TKGWV2 / KING-robust / READv2 kinship coefficients vs degree-threshold bands.
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
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=HERE / "input_values.tsv")
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

with args.input.open(newline="", encoding="utf-8") as handle:
    values = {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}

focal = "#c1121f"
methods=["TKGWV2","KING-robust","READv2"]
coeffs=[float(values["TKGWV2_half_relatedness"]),
        float(values["KING_robust_kinship"]),
        round(float(values["READv2_normalized_kinship"]), 3)]
snps=["4.5M SNPs","5.3M SNPs","169k SNPs"]
king_kin=float(values["KING_robust_kinship"])
ibs0_obs=float(values["IBS0_observed_per_site"])
ibs0_fs=float(values["IBS0_full_sibling_expectation"])
ibs0_unrel=float(values["IBS0_unrelated_expectation"])
ref_pts={"MZ twin / duplicate":(0.50,0.0),"Parent-offspring":(0.25,0.0),
         "Full siblings":(0.25,ibs0_fs),"2nd degree":(0.125,ibs0_fs*2),"Unrelated":(0.0,ibs0_unrel)}
fig,(axA,axB)=plt.subplots(1,2,figsize=(7.6,3.7),gridspec_kw={"width_ratios":[1.0,1.15]})
for lo,hi,col in [(0.177,0.354,"#c7e9c0"),(0.0884,0.177,"#f0f0f0"),(0.0,0.0884,"#ffffff")]:
    axA.axhspan(lo,hi,color=col,zorder=0)
axA.axhline(0.177,color="#31a354",lw=0.9,ls="--"); axA.axhline(0.0884,color="#bdbdbd",lw=0.6,ls=":")
axA.text(2.54,0.315,"1st-degree\nband",ha="right",va="center",fontsize=6,color="#31a354")
x=np.arange(3); axA.scatter(x,coeffs,s=70,color=focal,zorder=5,edgecolor="black",linewidth=0.5)
for xi,c,s in zip(x,coeffs,snps):
    axA.annotate(f"{c:.3f}",(xi,c),textcoords="offset points",xytext=(0,8),ha="center",fontsize=6.5,fontweight="bold")
    axA.annotate(s,(xi,c),textcoords="offset points",xytext=(0,-13),ha="center",fontsize=5.5,color="#555555")
axA.annotate("type: parent-offspring",(2,0.231),textcoords="offset points",xytext=(0,-23),ha="center",fontsize=5.5,color=focal,fontstyle="italic")
axA.set_xticks(x); axA.set_xticklabels(methods, fontsize=7.5); axA.set_ylim(0,0.40); axA.set_xlim(-0.6,2.6)
axA.set_ylabel("Kinship coefficient",labelpad=10); axA.set_title("Three methods agree: first-degree",loc="left",fontsize=8)
axA.annotate("a",xy=(0,1),xycoords="axes fraction",xytext=(-42,10),textcoords="offset points",fontsize=11,fontweight="bold",va="bottom",annotation_clip=False)
for name,(kx,ky) in ref_pts.items():
    axB.scatter(kx,ky,s=45,facecolor="white",edgecolor="#888888",linewidth=0.9,zorder=3)
axB.scatter(king_kin,ibs0_obs,s=95,color=focal,edgecolor="black",linewidth=0.6,zorder=6)
axB.annotate("ORT15 x ORT16\n(0.238, 0.0034)",(king_kin,ibs0_obs),textcoords="offset points",xytext=(-6,10),ha="right",fontsize=6,color=focal,fontweight="bold")
lab_off={"MZ twin / duplicate":(0,9,"center","bottom"),"Parent-offspring":(6,6,"left","center"),"Full siblings":(8,0,"left","center"),"2nd degree":(8,0,"left","center"),"Unrelated":(6,0,"left","center")}
for name,(kx,ky) in ref_pts.items():
    dx,dy,ha,va=lab_off[name]; axB.annotate(name,(kx,ky),textcoords="offset points",xytext=(dx,dy),ha=ha,va=va,fontsize=6,color="#555555")
axB.annotate("",xy=(0.25,ibs0_fs-0.001),xytext=(0.25,0.001),arrowprops=dict(arrowstyle="<->",color="#bbbbbb",lw=0.8))
axB.text(0.263,ibs0_fs/2,"PO vs FS\nseparate here",fontsize=5.5,color="#999999",va="center")
axB.set_xlim(-0.03,0.62); axB.set_ylim(-0.004,0.082)
axB.set_xlabel("KING kinship coefficient"); axB.set_ylabel("IBS0 rate (per site)",labelpad=8)
axB.set_title("IBS0 places the pair at parent-offspring",loc="left",fontsize=8)
axB.annotate("b",xy=(0,1),xycoords="axes fraction",xytext=(-42,10),textcoords="offset points",fontsize=11,fontweight="bold",va="bottom",annotation_clip=False)
fig.tight_layout()
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, bbox_inches="tight")
print(f"wrote {args.output}")
