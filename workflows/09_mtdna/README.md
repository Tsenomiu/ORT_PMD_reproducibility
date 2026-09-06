# 09 Mitochondrial haplogroup assignment

This corrected workflow accepts one canonical, collapsed-only, deduplicated
full-UDG BAM per individual. The BAM combines the independently deduplicated
collapsed reads from full-UDG libraries a, b and c. In the retained analysis this
is the `.1` state.

The `.2` state contains collapsed plus collapsed-truncated reads, and `.3` contains
those two classes plus paired reads. They are overlapping alternative read-set
states, not independent inputs. Never merge `.1`, `.2` and `.3`; doing so counts
collapsed records up to three times and collapsed-truncated records up to twice.
`validate_manifest.py` enforces one canonical collapsed-only row per sample and
rejects `.2`/`.3` filenames and non-canonical state labels. During execution it
also verifies the BAM against the known per-sample SHA-256, so renaming a nested
or otherwise different BAM cannot bypass the input guard.

## Calling workflow

`run_haplogroup.sh` extracts only the `MT` contig from the competitive hs37d5
alignment, calls haploid variants, reheaders the sample and runs HaploGrep 3. The
executed bcftools v1.16 call is:

```text
bcftools mpileup -r MT -f hs37d5.fa -B -q 30 -Q 30 MT.bam -Ou |
bcftools call --ploidy 1 -m
```

`-B` disables BAQ. Mapping quality and base quality are each required to be at
least 30. The bcftools v1.16 default maximum input depth of 250 reads per file is
retained because `-d` is not specified. No global minimum-depth,
allele-fraction or strand-balance filter is applied.

Reads were aligned competitively with BWA v0.7.17 to hs37d5, which contains
nuclear, decoy and rCRS mitochondrial sequence. Haplogroup classification uses
HaploGrep 3 v3.3.2 with `phylotree-fu-rcrs@1.3`. HaploGrep scores are software
outputs interpreted together with the underlying call set; this workflow does
not impose a universal score threshold.

## Inputs and execution

Edit a copy of `full_udg_bam_manifest.example.tsv`, retaining its four-column
schema and the canonical state/scope values. Then run each sample independently:

```bash
export ORT_CONFIG="$PWD/workflows/config.sh"
workflows/09_mtdna/run_haplogroup.sh full_udg_bams.tsv ORT15 work/mtdna
workflows/09_mtdna/run_haplogroup.sh full_udg_bams.tsv ORT16 work/mtdna
```

The private BAMs, VCFs and HaploGrep runtime are not distributed. Path-neutral
input and output checksums are provided in `corrected_input_checksums.tsv` and
`corrected_output_checksums.tsv`; tool and reference provenance are recorded in
`environment.tsv`.

## Pair comparison and target-site audit

`summarize_target_sites.py` reports aggregate support at selected sites under the
MAPQ/base-quality filters without emitting read names or private BAM paths.
`summarize_calls.py` compares called SNPs by chromosome, position, REF and ALT,
excludes the artificial `MT:3106 CN>C` rCRS-spacer record from biological SNP
counts, and keeps “no qualifying coverage” distinct from a reference genotype.

The corrected linear-rCRS call sets contain 36 ORT15 SNPs and 37 ORT16 SNPs, with
36 exact shared calls. The sole additional ORT16 caller-emitted SNP is C7028T,
supported by one forward read; it is not treated as a confirmed pairwise
difference. D4o1 assignments and HaploGrep scores remain unchanged.

The 36/37 values are caller-emitted counts from the stated linear-rCRS workflow,
not complete breakpoint-independent mitogenome totals. A shifted/circular
reference sensitivity is appropriate if complete whole-mitogenome call counts
are required.
