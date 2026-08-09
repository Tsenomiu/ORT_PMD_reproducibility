# Supplementary Figure S3 — petrous micro-CT

The exporter documents the deterministic display procedure. Native VOX volumes and
slice coordinates are not distributed, so this figure is not run by
`reproduce_figures.sh`.

```bash
python figures/supplementary_03_ct/export_ct_figure.py \
  ORT15.VOX X15 Y15 Z15 ORT16.VOX X16 Y16 Z16 > figure.png
```

`SCAN_METADATA.md` records the acquisition and display settings. The exporter
selects recorded slices, applies the shared display window, averages five adjacent
slices, and adds the scale bars; it does not alter the depicted anatomy.
