#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
# Publication wrapper: fetch the pinned GPL-2.0 TKGWV2 source and apply the
# exact compatibility patch recovered from the ORT analysis environment.
set -euo pipefail

readonly UPSTREAM_URL="https://github.com/danimfernandes/tkgwv2.git"
readonly UPSTREAM_COMMIT="c8638d47d3143b82ec969259df66e16f63b2b0ae"

usage() {
  echo "Usage: $0 OUTPUT_DIRECTORY" >&2
}

if [[ $# -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  usage
  exit 0
fi

if [[ $# -ne 1 ]]; then
  usage
  exit 2
fi

output_dir=$1
[[ ! -e "$output_dir" ]] || {
  echo "Refusing to overwrite existing path: $output_dir" >&2
  exit 3
}

for command_name in git; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "Required command not found: $command_name" >&2
    exit 4
  }
done

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    echo "Neither sha256sum nor shasum is available" >&2
    exit 4
  fi
}

check_hash() {
  local expected=$1
  local file=$2
  local observed
  observed=$(sha256_file "$file")
  [[ "$observed" == "$expected" ]] || {
    echo "SHA-256 mismatch: $file" >&2
    echo "Expected: $expected" >&2
    echo "Observed: $observed" >&2
    exit 5
  }
}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
patch_file="$script_dir/TKGWV2_compatibility.patch"
check_hash c562fd0d98dfa13cd14262a98a0b7b85e44be0b3c3796017ff4c543259cee7b8 "$patch_file"

git clone --quiet --no-checkout -- "$UPSTREAM_URL" "$output_dir"
git -C "$output_dir" checkout --quiet --detach "$UPSTREAM_COMMIT"

check_hash bb21bc649a79683367cc3fbc01e434fb64bafb412620de6d6ef5bc27845cc5cf "$output_dir/TKGWV2.py"
check_hash f75c5e5b9fd78897ccd5fec1f07b78a80bda106d1d2f64e875d3fd9023d76c58 "$output_dir/scripts/pileup2ped.py"
check_hash 9c215f7be6521f3c7bd524e9cb5b52bd1939e4562e2375124ac1a07a275e3639 "$output_dir/scripts/bam2plink.R"
check_hash 93e891c00c05273ad76a93699bbe91d8b0a68755b39835de74315ee38289d307 "$output_dir/scripts/plink2tkrelated.R"

git -C "$output_dir" apply "$patch_file"

check_hash 9187f573d4f6fe99af06b68295e1af48fd4fa9f56d742b6c5536b056ff73388b "$output_dir/scripts/bam2plink.R"
check_hash bb9a20e4ff114a051a111716d913fc6e42fe8e2dbf3528aa3bf14ff513b30ee3 "$output_dir/scripts/plink2tkrelated.R"

echo "Prepared TKGWV2 v1.0b at $UPSTREAM_COMMIT with the recovered compatibility patch."
