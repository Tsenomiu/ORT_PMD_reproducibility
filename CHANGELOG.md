# Changelog

## Unreleased

- Remove server-specific identifiers from the privacy checker, scan the checker
  itself, and report matching filenames without disclosing their contents.

## v1.2.1

- Shorten the README and consolidate release reports into one validation guide.
  Historical reports remain available in their original Git tags.
- Remove the duplicate limitations document.
- Align supplementary figure directories, output names and validation labels
  with the manuscript: S1–S4 and S7–S8 are rebuilt; S5–S6 require external inputs.
- Keep a checksum-pinned TKGWV2 expected table separate from the figure source
  and reproduced results, with regression tests that reject altered values.
- No scientific results or analysis parameters change.

## v1.2.0

- Add the Reduce18F/Round45 DeltaALT uncertainty analysis: 26 LOCO state
  intervals, 90 primary paired contrasts, 24 common-callable sensitivity
  contrasts, path-neutral provenance scripts and Figure 2 confidence intervals.
- Validate correction-versus-uncorrected effects under paired LOCO, 5-Mb and
  10-Mb genomic-block jackknives and document the external per-site TSV boundary.

- Correct the mitochondrial workflow to use one canonical collapsed-only,
  deduplicated full-UDG BAM per sample and reject overlapping `.2`/`.3` nested
  read-set states.
- Add checksum-enforced manifest validation, aggregate target-site support, exact-call
  comparison code, synthetic regression tests and checksummed corrected
  provenance.
- Replace the defective 37/37 mitochondrial summary with 36 ORT15 and 37 ORT16
  caller-emitted SNPs, of which 36 match exactly by position, REF and ALT.
- Retain D4o1 and HaploGrep scores 0.8911/0.9107; remove the unsupported ORT15
  A4215G claim and classify one-read ORT16 C7028T as unconfirmed.
- Record `MT:3106 CN>C` in both samples as the non-biological artificial rCRS
  N-spacer deletion and keep it outside biological SNP/discordance counts.
- Retire the obsolete Figure 1 SVG and rendering scripts. Figure 1 is an
  author-prepared PowerPoint schematic and is no longer claimed as a
  repository-generated figure.
- Preserve every v1.1.x historical release record without revision.

## v1.1.1

- Add an explicit `ci-portable` PCA raster-validation mode for GitHub Actions.
- Retain `canonical-strict` as the default, including all 20 PCA checks and all
  three exact canonical raster SHA-256 comparisons.
- Document the renderer-dependent anti-aliasing difference between canonical
  Poppler 25.06.0 / TeX Live 2024 output and Poppler 24.02.0 in the failing
  v1.1.0 GitHub Actions runs.
- No scientific data, analysis results, parameters, source values, figure
  content, tables, or conclusions changed.

## v1.1.0

- Add the recovered TKGWV2 v1.0b compatibility patch, a path-neutral
  processed-input runner, input and environment manifests, and exact validation
  of all 16 reported comparisons.
- Add READv2 downstream commands, normalization data, and validators for the
  primary, autosome-only, and transversion-only calculations, with the missing
  primary BAM-to-genotype command and random state disclosed.
- Add the recovered TEST/RC PCA preparation, projection, metric, and plotting
  records in publication-safe form, together with 20 validation checks and
  three exact raster comparisons.
- Add the recovered original Table 1 generator, its generated output, and the
  16-row/80-cell provenance check.
- Add main-output source-data, workflow, validation, file, and exclusion
  manifests for the proposed release.

## v1.0.1

- Simplified the reproducibility package while preserving the reported results
  and supported figure regeneration.
- Added compact KING/IBS0 and READv2 summary evidence and improved portability
  checks.
