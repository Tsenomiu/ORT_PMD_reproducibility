# PCA-derived plotting inputs

`query_coordinates.tsv` contains the 48 ORT projections used in Figure 7. Axis
variance percentages are calculated from
`data/summary/pca/pca_aadr_matched48.eval`.

`reference_aggregates.npz` is a disclosure-controlled plotting derivative of
the licensed AADR-v62 reference projection. For each continent and panel it
stores the plotting grid, a boolean 92%-coverage KDE mask, and fixed-grid bin
centres/counts only for cells with at least three contributors. It also stores
the Yakut PC1/PC2 mean and sample count and the fixed zoom window. It contains
no reference-individual IDs, labels, or coordinates and cannot reproduce the
upstream smartpca analysis.
