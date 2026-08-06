# 04 Reference and alternative allele counts

`run_alt_counts.sh` counts REF and ALT reads at 1000 Genomes biallelic sites with
bcftools, MAPQ at least 30, base quality at least 20, indels excluded, and a depth
cap of 10,000.

`summarize_alt_fraction.py` uses two explicit denominators:

- mean ALT fraction, `ALT/(REF+ALT)`, across sites with analysed depth at least 3;
- retained-site counts across every covered site, relative to the raw treatment
  for the same individual and library.

C-to-T and G-to-A sites form the damage-transition group. Counts for every
treatment must come from that treatment's own pileup.
