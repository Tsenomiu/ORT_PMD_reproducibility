# 04 Alternate-read fractions and site retention

`run_alt_counts.sh` counts REF and ALT reads at 1000 Genomes biallelic sites with
bcftools, MAPQ at least 30, base quality at least 20, indels excluded, and a depth
cap of 10,000.

`summarize_alt_fraction.py` defines analysed depth as `AD_REF + AD_ALT` and uses
two explicit thresholds:

- mean ALT fraction, `AD_ALT/(AD_REF+AD_ALT)`, across sites with analysed depth
  at least 3;
- retained-site counts across sites with analysed depth greater than 0, relative
  to the uncorrected library of the same individual and library type.

C-to-T and G-to-A sites form the damage-transition group. Every treatment is
summarised with the same rules and its own callable site set; the workflow does
not impose a fixed common-site intersection. Counts for every treatment must come
from that treatment's own pileup.
