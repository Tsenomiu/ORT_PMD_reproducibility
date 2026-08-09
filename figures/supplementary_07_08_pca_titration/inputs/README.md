# PCA-titration-derived plotting inputs

`query_coordinates.tsv` contains the 1,694 published ORT titration projections
used in Supplementary Figures S7 and S8.

`reference_aggregates.npz` is a disclosure-controlled plotting derivative of
the licensed AADR-v62 reference projection. For each individual/continent
window it stores a plotting grid, a boolean 92%-coverage KDE mask, fixed-grid
bin centres/counts only for cells with at least three contributors, and a
non-identifying visibility flag. It also stores the Yakut PC1/PC2 mean and
sample count and the two fixed plotting windows. It contains no
reference-individual IDs, labels, or coordinates and cannot reproduce the
upstream smartpca analysis.
