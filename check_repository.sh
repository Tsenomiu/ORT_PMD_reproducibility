#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
PYTHON=${PYTHON:-python3}

PASS=0
SKIP=0

pass() {
  printf 'PASS  %s\n' "$1"
  PASS=$((PASS + 1))
}

skip() {
  printf 'SKIP  %s\n' "$1"
  SKIP=$((SKIP + 1))
}

fail() {
  printf 'FAIL  %s\n' "$1" >&2
  exit 1
}

TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/ort-pmd-validate.XXXXXX")
trap 'rm -rf "$TMP_ROOT"' EXIT

# Shell syntax
while IFS= read -r -d '' script; do
  bash -n "$script" || fail "shell syntax: $script"
done < <(find . -type f -name '*.sh' -not -path './_build/*' -print0)
pass 'shell syntax'

# Python syntax without writing caches into the release tree.
if command -v "$PYTHON" >/dev/null 2>&1; then
  while IFS= read -r -d '' script; do
    PYTHONPYCACHEPREFIX="$TMP_ROOT/pycache" "$PYTHON" -m py_compile "$script" || \
      fail "Python syntax: $script"
  done < <(find . -type f -name '*.py' -not -path './_build/*' -print0)
  pass 'Python syntax'
else
  skip "Python syntax ($PYTHON unavailable)"
fi

# Excluded human-genomic, reference, raw-image, and runtime files. Capture the first
# result instead of piping `find -quit` into `grep -q`: under `pipefail`, grep can
# close the pipe early and turn a real hit into a false negative via SIGPIPE.
BAD_FILE=$(find . -type f \( \
    -name '*.bam' -o -name '*.bai' -o -name '*.cram' -o -name '*.crai' -o \
    -name '*.sam' -o -name '*.vcf' -o -name '*.vcf.gz' -o -name '*.vcf.bgz' -o \
    -name '*.bcf' -o -name '*.tbi' -o -name '*.csi' -o -name '*.h5' -o \
    -name '*.hdf5' -o -name '*.bed' -o -name '*.bim' -o -name '*.fam' -o \
    -name '*.pgen' -o -name '*.pvar' -o -name '*.psam' -o -name '*.geno' -o \
    -name '*.ind' -o -name '*.snp' -o -name '*.fastq' -o -name '*.fastq.gz' -o \
    -name '*.fq' -o -name '*.fq.gz' -o -name '*.fa' -o -name '*.fasta' -o \
    -iname '*.vox' -o -iname '*.dcm' -o -iname '*.dicom' -o \
    -name '*.pyc' -o -name '.DS_Store' \
  \) -not -path './_build/*' -print -quit)
if [[ -n "$BAD_FILE" ]]; then
  find . -type f \( \
    -name '*.bam' -o -name '*.bai' -o -name '*.cram' -o -name '*.crai' -o \
    -name '*.sam' -o -name '*.vcf' -o -name '*.vcf.gz' -o -name '*.vcf.bgz' -o \
    -name '*.bcf' -o -name '*.tbi' -o -name '*.csi' -o -name '*.h5' -o \
    -name '*.hdf5' -o -name '*.bed' -o -name '*.bim' -o -name '*.fam' -o \
    -name '*.pgen' -o -name '*.pvar' -o -name '*.psam' -o -name '*.geno' -o \
    -name '*.ind' -o -name '*.snp' -o -name '*.fastq' -o -name '*.fastq.gz' -o \
    -name '*.fq' -o -name '*.fq.gz' -o -name '*.fa' -o -name '*.fasta' -o \
    -iname '*.vox' -o -iname '*.dcm' -o -iname '*.dicom' -o \
    -name '*.pyc' -o -name '.DS_Store' \
  \) -not -path './_build/*' -print >&2
  fail 'excluded binary/genomic/runtime files found'
fi
pass 'excluded-file scan'

# Private infrastructure. The checker excludes itself because it contains this
# detection expression; all other text files must be portable.
PRIVATE_PATTERN='/media/test/|/home/(rc|test)/|/Users/[^/]+/ORT/|133[.]28[.]62[.]238|ssh[[:space:]]+(test|rc)([[:space:]]|$)'
if command -v rg >/dev/null 2>&1; then
  if rg -n --hidden "$PRIVATE_PATTERN" . \
      --glob '!check_repository.sh' --glob '!_build/**'; then
    fail 'private host or absolute working path found'
  fi
else
  if grep -RInE "$PRIVATE_PATTERN" . \
      --exclude='check_repository.sh' \
      --exclude-dir='_build' --exclude='*.pdf'; then
    fail 'private host or absolute working path found'
  fi
fi
pass 'private-path scan'

# High-confidence credential forms. Generic words such as "token" are documented
# elsewhere and are not themselves evidence of a leaked secret.
SECRET_PATTERN='BEGIN (RSA |OPENSSH |EC |DSA )?PRIVATE KEY|AKIA[0-9A-Z]{16}|ASIA[0-9A-Z]{16}|xox[baprs]-[A-Za-z0-9-]+'
if command -v rg >/dev/null 2>&1; then
  if rg -n --hidden "$SECRET_PATTERN" . \
      --glob '!check_repository.sh' --glob '!_build/**'; then
    fail 'credential-like material found'
  fi
else
  if grep -RInE "$SECRET_PATTERN" . \
      --exclude='check_repository.sh' \
      --exclude-dir='_build' --exclude='*.pdf'; then
    fail 'credential-like material found'
  fi
fi
pass 'credential scan'

# Minimal CFF structural check; use a YAML parser as an additional check when present.
grep -q '^cff-version: 1[.]2[.]0$' CITATION.cff || fail 'CITATION.cff version'
grep -q '^title:' CITATION.cff || fail 'CITATION.cff title'
grep -q '^authors:' CITATION.cff || fail 'CITATION.cff authors'
if command -v "$PYTHON" >/dev/null 2>&1 && "$PYTHON" -c 'import yaml' >/dev/null 2>&1; then
  "$PYTHON" -c 'import yaml; yaml.safe_load(open("CITATION.cff", encoding="utf-8"))' || \
    fail 'CITATION.cff YAML parse'
  pass 'CITATION.cff structure and YAML parse'
else
  pass 'CITATION.cff required fields (PyYAML unavailable)'
fi

for required in README.md LICENSE CITATION.cff requirements.txt workflows figures data/summary docs; do
  [[ -e "$required" ]] || fail "required path missing: $required"
done
pass 'required repository structure'

# Verify the release file manifest after all structural and content checks. The
# manifest excludes itself and runtime products under _build/.
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c MANIFEST.sha256 >/dev/null || fail 'checksum manifest'
  pass 'checksum manifest'
elif command -v shasum >/dev/null 2>&1; then
  shasum -a 256 -c MANIFEST.sha256 >/dev/null || fail 'checksum manifest'
  pass 'checksum manifest'
else
  skip 'checksum manifest (no SHA-256 utility available)'
fi

printf '\nRepository checks complete: %d PASS, %d SKIP, 0 FAIL\n' "$PASS" "$SKIP"
