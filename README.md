# geneset-modulescoring

A function to calculate gene module activity based on the approach described by [Tirosh _et al._](10.1126/science.aad0501). This is a repurposing of the _AddModuleScore()_ function from [Seurat](https://satijalab.org/seurat/index.html), adapted for bulk RNA sequencing data processed with [DESeq2](http://dx.doi.org/10.1186/s13059-014-0550-8).

## Installation

Copy the contents of `geneset_modulescoring.R` into your R session or save it locally to call using `source()`. This is demonstrated in the tutorial notebooks.

## Description

The `AddGeneSetScore` function scores samples scored based on the average normalized expression of the genes within a gene set of interest. From this, the enrichment score of randomly selected control genes with similar expression levels are subtracted. These control gene sets are defined by first binning all genes into bins of aggregate expression level and then, for each gene in the gene set of interest, 100 genes from the same expression bin as that gene are randomly selected. In this way, the control gene sets have a comparable distribution to the gene set of interest and the control gene set is 100 times larger such that the score for the gene set of interest is analogous to averaging over 100 randomly-selected gene sets as the gene set of interest. The score is then set to range from zero, meaning no enrichment compared to random sets of genes with similar expression, to one, reflecting the highest average expression of all genes within the gene set of interest.

Input arguments are:

(_required_)
- `dds`: a DESeqDataSet
- `features`: a character vector of gene names

(_optional_)
- `pool`: list of features to check expression levels against, defaults to `rownames(x = counts(dds))`.
- `nbin`: number of bins of aggregate expression levels for all analyzed features
- `ctrl`: number of control features selected from the same bin per analyzed feature
- `name`: name for the expression program
- `seed`: set a random seed. If `NULL`, no seed is set.

If an input feature is not present in the rowdata of the DESeqDataSet, a warning message will appear and the module score will be calculated using the remaining features (cf. tutorial).

## Citation

When using this function, please cite: __reference coming soon__
