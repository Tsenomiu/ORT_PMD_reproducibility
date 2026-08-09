# Figure code for the associated article

`FIGURE_MANIFEST.tsv` maps the seven main and eight supplementary figures to their
deterministic generators and input status.

```bash
python -m pip install -r ../requirements-figures-lock.txt
bash figures/reproduce_figures.sh _build/figures
```

The script rebuilds all figures supported by included inputs. It reports explicit
skips for Supplementary Figure S2 (controlled archaeological coordinates and contour
source) and Supplementary Figure S3 (controlled VOX volumes). PCA backgrounds use
population-level aggregates; exact individual-level AADR backgrounds require
separately licensed inputs.

The main Figure 1 schematic is stored as editable SVG. Main Figures 2--7 and
Supplementary Figures S1, S4--S8 are generated from data or archived coordinates.
