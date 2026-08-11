# Curated READv2 processed data

These small derived tables support numerical verification without distributing
cemetery BAMs or complete genotype matrices.

- `analysis_summary.tsv`: generated publication summary of the three accepted
  ORT15–ORT16 analyses, including exact and three-decimal values.
- `focal_results.tsv`: processed source data containing the exact archived focal
  READv2 output fields.
- `sample_manifest.tsv`: processed source data listing all 27 primary input
  labels, the six `--mind 0.8` exclusions, the 21 retained labels, and all 22
  labels in the later autosome/transversion run.
- `primary_p0_anonymized.tsv`: all 210 primary pairwise raw and normalized P0
  values, numerically ranked with pair identities removed.
- `autosome_1240k_p0_anonymized.tsv`: all 230 later all-1240k/autosome pairwise
  P0 values, numerically ranked with pair identities removed.
- `transversion_only_p0_anonymized.tsv`: all 230 matching transversion-only
  pairwise P0 values, numerically ranked with pair identities removed.
- `validation_results.tsv`: generated output from
  `workflows/08_kinship/readv2/validate_reported_results.py`.

The anonymized tables retain the values and overlap counts needed to recompute
each within-cohort median, the ORT15–ORT16 normalized P0 and the kinship
coefficient. Unrelated cemetery pair identities are not duplicated here.

The complete primary 210-pair output already present under
`data/summary/kinship/` is not required by this validator. No corresponding full
pair-labelled later table is included.

See `workflows/08_kinship/readv2/README.md` for the exact software revision,
commands, normalization formula, rerun boundary and limitation.
