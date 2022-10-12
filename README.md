# SARS-CoV-2 mutational spectrum

## Overview
This repository analyzes the mutational spectrum of human SARS-CoV-2.

This analysis was originally based off starting to replicate analyses [Neher (2022)](https://www.biorxiv.org/content/10.1101/2022.08.22.504731v1.full), but then checking the assumption that underlying synonymous mutation spectra were conserved across clades.

Study by Jesse Bloom, Annabel Beichman, Richard Neher, Kelley Harris.

## Structure of repository and running the analysis
The analysis is entirely reproducible.

First build the [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml).
This requires you to install [conda](https://docs.conda.io/), and then run:

    conda env create -f environment.yml

That command will create a `conda` environment named `SARS2-mut-spectrum` which you can activate with:

    conda activate SARS2-mut-spectrum

Then run the [snakemake](https://snakemake.readthedocs.io/) pipeline in [Snakefile](Snakefile), which reads its configuration from [config.yaml](config.yaml) by running:

    snakemake -j <n_cpus> --use-conda

where `<n_cpus>` is the number of CPUs to use.
You need the `--use-conda` command because although almost all the pipeline is run by the `SARS2-mut-spectrum` `conda` environment in [environment.yml](environment.yml), the drawing of the tree using [ete3](http://etetoolkit.org/) uses the `conda` environment in [environment_ete3.yml](environment_ete3.yml).

Note that the pipeline uses Python scripts in [./scripts/](scripts) and Jupyter notebooks in [./notebooks/](notebooks).

The created files are placed in [./results/](results).
Only some of those results files are tracked in this repo (others are too large to track).

Interactive versions of some figures are rendered via GitHub Pages at [https://jbloomlab.github.io/SARS2-mut-spectrum/](https://jbloomlab.github.io/SARS2-mut-spectrum/)

## Overview of analysis

### Counting mutations
Mutation counts are extracted from the pre-built mutation-annotated tree that is made available for use with the [UShER](https://usher-wiki.readthedocs.io/) package.
This tree contains all public access SARS-CoV-2 sequences, with mutations annotated on branches.
For the specific version of the tree used here, see the [config.yaml](config.yaml) file.

We analyze mutations grouping sequences at the level of [Nextstrain clades](https://clades.nextstrain.org/), which are already annotated on the pre-built mutation-annotated tree.
For each Nextstrain clade, we use the clade founder genotype manually defined by [Neher (2022)](https://www.biorxiv.org/content/10.1101/2022.08.22.504731v1.full) and available at the URL indicated in [config.yaml](config.yaml).

We perform some crucial filtering to remove spurious mutations as can arise from bad sequencing, base calling to reference, etc:

 1. We ignore all mutations on any branches with high numbers of total mutations, or mutations to either the reference or the founder sequence of the clade in question. The exact settings for this filtering are in [config.yaml](config.yaml).

 2. We ignore any branches with multiple mutations at the same codon.

 3. If a switch is set in [config.yaml](config.yaml) (it currently is), we specify to exclude any mutations that are reversions from the clade founder to the reference, and also the reverse complement of these mutations. This is designed to remove missing bases called to reference, and also complements of those mutations induced by spurious nodes with such miscalls on downstream branches in the tree.

 4. We specify to exclude all mutations at error-prone or problematic sites as manually specified in [config.yaml](config.yaml).

 5. We ignore any clades with small numbers of sequence samples as indicated in [config.yaml](config.yaml) as these are expected to have too much noise.

 6. Note also that indels are ignored, as they are not captured in the mutation-annotated tree.

The above mutation counts both for all sequences for a clade, and for the sample subsets defined in [config.yaml](config.yaml) are stored in [results/mutation_counts/aggregated.csv](results/mutation_counts/aggregated.csv).
Note that mutations are annotated by the protein(s) they affect, if they are synonymous, if they are at 4-fold degenerate sites, etc.

### Analysis of 4-fold synonymous mutation spectrum / rates

We then analyze the mutation spectra and rate of mutations.
For this analysis, we only consider synonymous mutations at sites (third codon positions) that are four-fold degenerate in the founder sequence for each clade.
The file [results/clade_founder_nts/clade_founder_nts.csv](results/clade_founder_nts/clade_founder_nts.csv) specifies which these sites are for each clade.

Specifically, we determine the relative fraction of all 4-fold synonymous mutations that are each type of nucleotide change, and also the relative **rates** of the different types of mutations, which are just computed as the fraction of all mutations of that type normalized by the composition of the sequence at these 4-fold degenerate sites.
These rates for each clade with sufficient counts are written to the file [results/synonymous_mut_rates/rates_by_clade.csv](results/synonymous_mut_rates/rates_by_clade.csv).

We also perform analyses for subsets of sequences from different regions (as specified in [config.yaml](config.yaml)) as well as for the the genome partitioned into halves--these analyses are designed to check that results are not due to artifacts related to sequencing pipelines or hotspots in the genome.
For all of these analyses, we only include subsets/partitions with at least the minimum number of mutations indicated in [config.yaml](config.yaml).

Most of the analysis of the synonymous mutation spectrum is done by [notebooks/synonymous_mut_rates.ipynb](notebooks/synonymous_mut_rates.ipynb).
You can look at the HTML rendering of running that Jupyter notebook at [results/synonymous_mut_rates/synonymous_mut_rates.html](results/synonymous_mut_rates/synonymous_mut_rates.html): download and open that notebook to look at the interactive `altair` plots.

### Drawing of tree

We draw a tree of the clade founder sequences.
The tree itself is in [results/clade_founder_tree/clade_founders.treefile](results/clade_founder_tree/clade_founders.treefile), and a matrix of the distances between clade founders is in [results/clade_founder_tree/clade_founders.mldist](results/clade_founder_tree/clade_founders.mldist).

The notebook [draw_tree_w_mut_enrichments.ipynb](draw_tree_w_mut_enrichments.ipynb), which requires the `conda` environment in [environment_ete3.yml](environment_ete3.yml), draws the tree in [results/clade_founder_tree/tree_w_enrichments.png](results/clade_founder_tree/tree_w_enrichments.png), which has nodes plotted with change in relative mutation rate (absolute difference relative to reference clade, which is first listed).

### Mantel test
A Mantel test for correlation between mutation rates and phylogenetic signal is run using the Euclidean distances between the relative mutation rates in [results/synonymous_mut_rates/clade_rate_distances.csv](results/synonymous_mut_rates/clade_rate_distances.csv) and the phylogenetic distances in [results/clade_founder_tree/clade_founders.mldist](results/clade_founder_tree/clade_founders.mldist).

## Caveats of analysis
None of these are expected to seriously affect the accuracy of the current analysis, but they could become problematic if the same analysis is applied to substantially more diverged clades:

 - For computing the rates of mutations from 4-fold degenerate synonymous sites, we do not exclude sites / mutations to exclude from the computation of the sequence composition that is used to normalize the counts of different mutations to rates. This would become a problem if you start to specify a very large number of sites to exclude in [config.yaml](config.yaml), or if the sequences become highly diverged from the reference (as we exclude reversions to reference and their complement).
 
 - Multiple mutations in same codon on a branch are excluded from analysis. This is not expected to have much effect, as this is rare.

 - Indels are ignored.

 - The way that codon positions are assigned for identify synonymous sites will fail if there are non-in-frame indels.

 - Four-fold synonymous sites are identified in the clade founder, which could lead to mis-identification if seuqences in a clade become highly diverged from the founder.
