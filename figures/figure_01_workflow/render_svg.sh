#!/bin/sh
set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SOURCE="$HERE/figure_01_workflow.svg"
OUTPUT=${1:?Usage: render_svg.sh OUTPUT.pdf}
mkdir -p "$(dirname -- "$OUTPUT")"

if command -v inkscape >/dev/null 2>&1; then
  inkscape "$SOURCE" --export-type=pdf --export-filename="$OUTPUT"
elif command -v rsvg-convert >/dev/null 2>&1; then
  rsvg-convert -f pdf -o "$OUTPUT" "$SOURCE"
elif command -v cairosvg >/dev/null 2>&1; then
  cairosvg "$SOURCE" -o "$OUTPUT"
else
  printf '%s\n' 'SKIP: install Inkscape, librsvg (rsvg-convert), or CairoSVG.' >&2
  exit 77
fi
