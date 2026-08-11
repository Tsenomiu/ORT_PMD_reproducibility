# TKGWV2 manuscript workflow

This directory is the publication copy of the recovered TKGWV2 analysis used
for the ORT15–ORT16 comparisons. It is intentionally limited to small,
path-neutral code, metadata and derived results. It does not contain BAM,
FASTQ, VCF, PED/MAP or 1000 Genomes reference files.

## Recovery and validation status

The executed TKGWV2 source and the exact processed PED/MAP inputs were recovered
from author-controlled project storage. A fresh rerun from those deterministic
processed inputs reproduced **all 16 reported results exactly (16/16)**. In
particular, the full-UDG × full-UDG bamRefine-5 comparison reproduced:

- HRC: **0.2356**;
- used SNPs: **4,517,214**;
- counts0: **319,067**;
- counts4: **4,198,147**;
- relationship class: **first degree**.

Run the compact public validation with:

```bash
python3 workflows/08_kinship/tkgwv2/validation/validate_results.py
```

The expected terminal message is `TKGWV2 validation: 16/16 PASS`.

## File provenance

| Path | Classification | Purpose |
|---|---|---|
| `TKGWV2_compatibility.patch` | original recovered patch | Exact four-line diff between upstream v1.0b and the executed local source |
| `fetch_and_patch_tkgwv2.sh` | publication wrapper | Fetches the pinned upstream commit, verifies it and applies the recovered patch |
| `run_tkgwv2_from_processed_peds.sh` | reconstructed file | Path-neutral runner reconstructed from successful logs; not claimed as the historical top-level script |
| `dyads/*.tsv` | original recovered inputs | Exact dyad definitions used for the three comparison groups |
| `input_checksums.tsv` | generated provenance record | Sizes and SHA-256 checksums for the deterministic external inputs |
| `parameters.tsv` | recovered/generated documentation | Recorded parameters and deterministic-boundary notes |
| `environment.tsv` | recovered/generated documentation | Original recorded and fresh-validation environments |
| `source_data/tkgwv2_pair_results.tsv` | processed source data | The 16 values reported in the manuscript package |
| `validation/reproduced_results.tsv` | generated validation output | Normalized results from the fresh processed-input rerun |
| `validation/validation_results.tsv` | generated validation output | Cell-by-cell comparisons of reported and reproduced values |
| `validation/validate_results.py` | publication validation code | Rechecks all 16 comparisons locally |
| `validation/verify_input_checksums.py` | publication validation code | Checks external deterministic inputs against `input_checksums.tsv` |

The private recovery area retains the historical source copies and logs. They
are not duplicated here because they contain machine-specific paths and are not
needed to execute the public workflow.

## Upstream identity and licensing

- Software: TKGWV2 v1.0b (released July 2022)
- Upstream repository: <https://github.com/danimfernandes/tkgwv2>
- Pinned upstream commit: `c8638d47d3143b82ec969259df66e16f63b2b0ae`
- Upstream author recorded in the source: Daniel Fernandes
- Upstream licence: GPL-2.0

The recovered patch modifies GPL-licensed upstream source and is distributed
under GPL-2.0-only; see `GPL-2.0-only.txt` and `THIRD_PARTY_NOTICE.md`. The
publication wrappers, validation code and ORT-derived tables remain under the
repository's MIT licence. The upstream source tree itself is not vendored.

## What changed relative to upstream

`TKGWV2_compatibility.patch` is the exact recovered diff. It makes four changes:

1. treats the BAM suffix as a literal string and returns one scalar sample ID;
2. adds `--threads 32` to the BAM-conversion PLINK call;
3. replaces a large fixed-string `grep` intersection with a keyed AWK
   intersection; and
4. replaces a second large `grep` frequency lookup with keyed AWK.

The two AWK substitutions were required because the system `ugrep` 7.5
implementation exceeded its pattern-complexity limit. The PLINK thread option is
a performance change. The suffix fix prevents regular-expression interpretation
of the filename suffix.

## Rebuild the patched software

The following creates a new detached checkout and refuses to overwrite an
existing directory:

```bash
TKGWV2_DIR="${PROJECT_ROOT}/software/tkgwv2-v1.0b-ort"
bash workflows/08_kinship/tkgwv2/fetch_and_patch_tkgwv2.sh "${TKGWV2_DIR}"
```

The wrapper verifies the four relevant upstream files before applying the patch
and verifies the two modified files afterwards. Network access and Git are
required.

## External input boundary

The deterministic rerun begins from 16 individual PED/MAP pairs and the EAS
allele-frequency file listed in `input_checksums.tsv`. These files total more
than 13 GB and are not suitable for GitHub. Verify an obtained copy before use:

```bash
python3 workflows/08_kinship/tkgwv2/validation/verify_input_checksums.py \
  --ped-map-dir "${INPUT_DIR}/ped_map" \
  --frequency-file "${REFERENCE_DIR}/1000GP_EAS_allchr_sorted.frq"
```

The PED/MAP files were created on hs37d5/GRCh37 at 1000 Genomes Phase 3
biallelic SNP sites. BAM conversion used `samtools mpileup -Q 30 -q 30 -B` and
TKGWV2's `pileup2ped.py`. That converter uses unseeded `random.choice` for
pseudo-haploid sampling. Consequently, a new BAM-level conversion is not
expected to reproduce the historical PED files byte-for-byte; the checksummed
PED/MAP files are the deterministic boundary for exact reproduction.

The allele-frequency input contains Phase 3 East Asian (EAS) frequencies. It is
third-party reference data and must be obtained under the provider's terms.

## Rerun from verified processed inputs

Install Python 3, R with `data.table`, and PLINK 1.9, then run:

```bash
bash workflows/08_kinship/tkgwv2/run_tkgwv2_from_processed_peds.sh \
  --tkgwv2 "${TKGWV2_DIR}/TKGWV2.py" \
  --frequency-file "${REFERENCE_DIR}/1000GP_EAS_allchr_sorted.frq" \
  --ped-map-dir "${INPUT_DIR}/ped_map" \
  --output-dir "${OUTPUT_DIR}/tkgwv2"
```

The runner creates isolated `fu_fu`, `nu_nu` and `nu_fu` directories, links
the corresponding PED/MAP files, copies the recovered dyad definitions and
executes the following command in each directory:

```text
python3 TKGWV2.py plink2tkrelated --freqFile EAS_FREQUENCY_FILE --dyads dyads.txt
```

It refuses to overwrite an existing output directory and validates the new
results against the 16 reported rows when complete. TKGWV2 itself uses no
random operation in this processed-PED stage.
