# Figure code

This directory contains one plotting workflow for each figure that can be generated
from distributable inputs. The compact manifest maps article figure numbers to code
directories and reproducibility status.

```bash
python -m pip install -r ../requirements.txt
./reproduce_figures.sh ../_build/figures
```

Generated PDFs are written to the requested output directory and are not committed to
Git. The script rebuilds all figures supported by included inputs and reports clear
skips for figures that require controlled data.

- The archaeological site map requires controlled spatial inputs and is not included
  as a code directory.
- The micro-CT exporter is supplied for method transparency, but the original VOX
  volumes and slice coordinates are not distributed.
- Figure 9 and Supplementary Figures S5–S6 use population-level reference aggregates.
  Exact individual-level AADR backgrounds require separately licensed inputs.

Python dependencies are listed in the top-level `requirements.txt`. Figure 1 also
requires Inkscape, `rsvg-convert`, or CairoSVG for SVG-to-PDF conversion. Poppler is
optional and is used only to inspect generated PDFs.
