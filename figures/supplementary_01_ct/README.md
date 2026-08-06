# Supplementary Figure S1 — petrous micro-CT

The exporter documents the deterministic display procedure. Native VOX volumes and
slice coordinates are not distributed, so this figure is not run by
`reproduce_figures.sh`.

```bash
python figures/supplementary_01_ct/export_ct_figure.py \
  ORT15.VOX X15 Y15 Z15 ORT16.VOX X16 Y16 Z16 > figure.png
```

`SCAN_METADATA.md` records the acquisition and display settings. The exporter does
not perform generative enhancement, denoising, inpainting, or anatomical alteration.
