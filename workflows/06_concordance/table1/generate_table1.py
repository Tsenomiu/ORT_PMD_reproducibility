#!/usr/bin/env python3
"""
regenerate_imputation_table.py — rebuild tab:imputation_summary from the fixed-validator
GLIMPSE2 concordance run (concordance_fixedval_20260713 on TEST).

Inputs (harvested to JSON):
  impute_gcsv.json    : {ind:{treatment:(NRD%, best-guess r2 overall, dosage r2 overall)}}
                        parsed from conc_<ind>_<t>.error.spl.txt.gz GCsV row
                        (0-based cols f16=NRD, f17=best-guess r2, f18=dosage r2)
  rsquare_all16.json  : per-MAF-bin rows from conc_<ind>_<t>.rsquare.grp.txt.gz
                        cols: bin, n, meanMAF, c4=best-guess r2, c5=dosage r2

Table columns (settled design; referee required dosage and best-guess reported separately):
  Sample | Strategy | dosage r2 (MAF>=0.1%) | best-guess r2 (MAF>=0.1%) | NRD (%)

r2 aggregation: count-weighted mean over MAF bins 1..8 (MAF>=0.1%); bin 0 (MAF<0.1%,
near-monomorphic) excluded. NRD is the GCsV overall non-reference discordance.

COLUMN IDENTITY NOTE: dosage r2 >= best-guess r2 always (dosage retains posterior
uncertainty). In rsquare.grp that is c5 (dosage) > c4 (best-guess); the earlier
post_summary.sh header labelling col4=dosage was inverted and is corrected here.
"""
import json, csv, sys

def wmean(bins, col, lo=1):
    b=[x for x in bins if x['bin']>=lo]
    return sum(x['n']*x[col] for x in b)/sum(x['n'] for x in b)

def main(gcsv_json, rsq_json, out_csv):
    gcsv=json.load(open(gcsv_json)); rs=json.load(open(rsq_json))
    order=["raw","trim5","trim10","rescale5","rescale10","rescaled","bamrefine5","bamrefine10"]
    disp={"raw":"Raw (uncorrected)","trim5":"Trim-5","trim10":"Trim-10","rescale5":"Rescale-5",
          "rescale10":"Rescale-10","rescaled":"Rescale-12","bamrefine5":"bamRefine-5","bamrefine10":"bamRefine-10"}
    recs=[]
    for ind in ["ORT15","ORT16"]:
        for t in order:
            nrd,bg_g,ds_g=gcsv[ind][t]
            k=f"{ind}_{t}"
            recs.append(dict(sample=ind, treatment=t, strategy=disp[t],
                dosage_r2_maf01=round(wmean(rs[k],'c5'),4),
                bestguess_r2_maf01=round(wmean(rs[k],'c4'),4),
                NRD_pct=round(nrd,2)))
    with open(out_csv,"w",newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=["sample","treatment","strategy","dosage_r2_maf01","bestguess_r2_maf01","NRD_pct"])
        w.writeheader(); w.writerows(recs)
    print(f"wrote {out_csv}: {len(recs)} rows")

if __name__=="__main__":
    a=sys.argv
    main(a[1] if len(a)>1 else "impute_gcsv.json",
         a[2] if len(a)>2 else "rsquare_all16.json",
         a[3] if len(a)>3 else "tab_imputation_fixedval_values.csv")
