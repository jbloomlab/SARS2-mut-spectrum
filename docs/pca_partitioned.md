---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: pca_partitioned
---

Principal component analysis on mutation spectrum, both with and without partitioning the genome.
Each point represents a different clade.
Note that although the points for the partitioned clades are slightly differently positioned, the shift is consistent across clades indicating that the differences between clades are **not** driven by just one genome partition.
The chart is interactive:

  - mouse over points for details
  - click on specific clades in the legend to highlight their points
  - clck on the "exclude top" legend to show only points +/- exclusions of top mutations. 
  - shift-click to make multiple selections

Note that the PCA itself is **not** re-calculated even when some points are de-selected for display.
