# 03 Post-mortem damage correction

This directory provides parameterised commands for the three tested approaches:

- **trimBam:** trim 5 bases from both ends of full-UDG reads and 5 or 10 bases
  from non-UDG reads.
- **mapDamage2:** rescale with terminal modelling lengths of 5 or 10 bases, or
  the program default of 12 bases, using `--rescale --merge-libraries`.
- **bamRefine v0.2.1:** mask sites near read ends with thresholds of 5 or 10
  bases and a GRCh37 1000 Genomes biallelic-SNP list.

`build_bamrefine_sites.sh` creates the required six-column site file
(`ID CHROM 0 POS REF ALT`). Run damage profiling separately for every treatment;
do not derive corrected-treatment summaries by transforming raw counts.
