---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: pca_exclude_top
---

Principal component analysis on mutation spectrum, both with and without excluding the top most abundant mutations in each clade.
Each point represents a different clade.
Note that although the clades +/- the top mutations are slightly differently positioned, the shift is consistent across clades indicating that the differences between clades are **not** driven by the top mutations.
The chart is interactive:

  - mouse over points for details
  - click on specific clades in the legend to highlight their points
  - clck on the "exclude top" legend to show only points +/- exclusions of top mutations. 
  - shift-click to make multiple selections

Note that the PCA itself is **not** re-calculated even when some points are de-selected for display.
