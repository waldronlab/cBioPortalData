
# cBioPortalData

<!-- start badges here -->

[![Travis Build
Status](https://travis-ci.org/waldronlab/cBioPortalData.svg?branch=master)](https://travis-ci.org/waldronlab/cBioPortalData)
[![Build
status](https://ci.appveyor.com/api/projects/status/42kd6prni3o0q50b?svg=true)](https://ci.appveyor.com/project/waldronlab/cbioportaldata)
<!-- end badges here -->

## cBioPortal data and MultiAssayExperiment

### Overview

This project aims to import all cBioPortal datasets as
[MultiAssayExperiment](http://bioconductor.org/packages/MultiAssayExperiment/)
objects in Bioconductor. It offers some advantages over the CDGS-R
package:

1.  The MultiAssayExperiment class explicitly links all assays to the
    patient clinical/pathological data
2.  The MultiAssayExperiment class provides a [flexible
    API](https://github.com/waldronlab/MultiAssayExperiment/wiki/MultiAssayExperiment-API)
    including harmonized subsetting and reshaping to convenient wide and
    long formats.
3.  It provides complete datasets, not just for subsets of genes
4.  It provides automatic local caching, thanks to BiocFileCache.

## MultiAssayExperiment Cheatsheet

<a href="https://github.com/waldronlab/cheatsheets/blob/master/MultiAssayExperiment_QuickRef.pdf"><img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/master/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/></a>

## Minor note

cBioPortalData is a work in progress. Due to variation in data
formatting, the package does not import all 268 (as of Dec 2019)
datasets. Currently, it successfully imports 93% of 200 randomly sampled
datasets. Please file an issue to request prioritization of any of the
remaining datasets.

## Quick Start

### Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("cBioPortalData", quietly = TRUE))
    BiocManager::install("cBioPortalData")

library(cBioPortalData)
```

### API Service

Flexible and granular access to cBioPortal data from
`cbioportal.org/api`. This option is best used with a particular gene
panel of interest. It allows users to download sections of the data with
molecular profile and gene panel combinations within a study.

``` r
cbio <- cBioPortal()
gbm <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "gbm_tcga",
    genePanelId = "IMPACT341",
    molecularProfileIds = c("gbm_tcga_rppa", "gbm_tcga_mrna")
)
```

``` r
gbm
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] gbm_tcga_rppa: SummarizedExperiment with 67 rows and 244 columns
#>  [2] gbm_tcga_mrna: SummarizedExperiment with 334 rows and 401 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save all data to files
```

### Packaged Data Service

This function will download a dataset from the `cbioportal.org/datasets`
website as a packaged tarball and serve it to users as a
`MultiAssayExperiment` object. This option is good for users who are
interested in obtaining all the data for a particular study.

``` r
laml <- cBioDataPack("laml_tcga")
```

``` r
laml
#> A MultiAssayExperiment object of 11 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 11:
#>  [1] CNA: SummarizedExperiment with 24776 rows and 191 columns
#>  [2] RNA_Seq_expression_median: SummarizedExperiment with 19720 rows and 179 columns
#>  [3] RNA_Seq_mRNA_median_Zscores: SummarizedExperiment with 19720 rows and 179 columns
#>  [4] RNA_Seq_v2_expression_median: SummarizedExperiment with 20531 rows and 173 columns
#>  [5] RNA_Seq_v2_mRNA_median_Zscores: SummarizedExperiment with 20531 rows and 173 columns
#>  [6] cna_hg19.seg: RaggedExperiment with 13571 rows and 191 columns
#>  [7] linear_CNA: SummarizedExperiment with 24776 rows and 191 columns
#>  [8] methylation_hm27: SummarizedExperiment with 10919 rows and 194 columns
#>  [9] methylation_hm450: SummarizedExperiment with 10919 rows and 194 columns
#>  [10] mutations_extended: RaggedExperiment with 2584 rows and 197 columns
#>  [11] mutations_mskcc: RaggedExperiment with 2584 rows and 197 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save all data to files
```
