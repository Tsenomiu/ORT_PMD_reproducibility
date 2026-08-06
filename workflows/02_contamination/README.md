# 02 Contamination estimation

`run_xchr_contamination.sh` runs the ANGSD v0.933 male-X analysis for ORT15 over
X:5,000,000–154,900,000 with MAPQ 30, base quality 20, and the HapMapChrX panel.

`run_mtdna_contamination.sh` runs the standalone mtCont component of schmutzi
v1.5.7 for the non-UDG libraries. It uses separate 5-prime and 3-prime damage
profiles, endoCaller likelihoods, and a single-contaminant-haplotype assumption.
This is a standalone mtCont analysis, not the complete iterative schmutzi workflow.

Input BAMs, likelihood matrices, HapMapChrX data, and the schmutzi Eurasian
frequency panel must be supplied locally.
