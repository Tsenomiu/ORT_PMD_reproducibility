# Figure 1 — workflow

This non-data schematic is stored as editable SVG and is synchronized with the
associated article. It should not be described as a data plot.

The editable SVG is rendered with the bundled Matplotlib script so that Figure 1
uses the same embedded DejaVu Sans typeface as Figures 2--7. The shell wrapper uses
the Python renderer when Matplotlib is available and otherwise falls back to a
standard SVG converter:

```bash
figures/figure_01_workflow/render_svg.sh _build/figures/figure_01_workflow.pdf
```
