
<!-- README.md is generated from README.Rmd. Please edit that file -->
timecourseRnaseq
================

<!-- badges: start -->
<!-- badges: end -->
The timecourseRnaseq package currently has one primary function:    
   
**check_enrichment()**    
   
which performs overrepresentation analyses (using hypergeometric test) and returns the enrichment results as a table. Additionally, the function also plots the enrichment results by default.


Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("biplabendu/timecourseRnaseq")
```

Usage
-------

> ?check_enrichment 

``` r
library(timecourseRnaseq)
some.beau.genes <- c("BBA_00100", "BBA_01000", "BBA_10000", "BBA_10001", "BBA_10002", "BBA_10003", "BBA_10004")

## GO enrichments + plot
check_enrichment(some.beau.genes, org="beau")

## pfam enrichments (no plot)
check_enrichment(some.beau.genes, org="beau", what="pfams", plot=F)

## pfam enrichments + plot + long-formatted enrichment results
check_enrichment(some.beau.genes, org="beau", what="pfams", plot=T, expand=T)
   # please note, expand=T works only if the annotation file contains a column 
   # named "gene_desc" that contains descriptions of genes

## specify your own annotation file (?check_enrichment)
some.speciesXX.genes <- c("genes", "of", "interest", "from", "speciesXX") # geneset
path <- "path/to/speciesXX/annotation/file.csv" 
check_enrichment(geneset=some.speciesXX.genes, path_to_annot = path, what="GOs") # double-check if sep needs to be specified

```
    
    
The annotation files for the four organisms (cflo_annots, beau_annots, ophio_cflo_annots, ophio_kim_annots) are provided.    
To see an annotation file, just type the following:

```
## org = "beau"
head(beau_annots)

## org = "cflo"
head(cflo_annots)

```
