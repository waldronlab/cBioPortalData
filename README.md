
# cBioPortalData

<!-- start badges here -->

[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/cBioPortalData.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/cBioPortalData)
[![Platforms](http://www.bioconductor.org/shields/availability/release/cBioPortalData.svg)](https://www.bioconductor.org/packages/release/bioc/html/cBioPortalData.html#archives)
[![Travis Build
Status](https://travis-ci.org/waldronlab/cBioPortalData.svg?branch=devel)](https://travis-ci.org/waldronlab/cBioPortalData)
[![Build
status](https://ci.appveyor.com/api/projects/status/42kd6prni3o0q50b?svg=true)](https://ci.appveyor.com/project/waldronlab/cbioportaldata)
[![Downloads](https://www.bioconductor.org/shields/downloads/release/cBioPortalData.svg)](https://bioconductor.org/packages/stats/bioc/cBioPortalData/)
<!-- end badges here -->

## cBioPortal data and MultiAssayExperiment

### Overview

The `cBioPortalData` R package aims to import cBioPortal datasets as
[MultiAssayExperiment](http://bioconductor.org/packages/MultiAssayExperiment/)
objects into Bioconductor. Some of the features of the package include:

1.  The use of the `MultiAssayExperiment` integrative container for
    coordinating and representing the data.
2.  The data container explicitly links all assays to the patient
    clinical/pathological data.
3.  With a [flexible
    API](https://github.com/waldronlab/MultiAssayExperiment/wiki/MultiAssayExperiment-API),
    `MultiAssayExperiment` provides harmonized subsetting and reshaping
    into convenient wide and long formats.
4.  The package provides datasets from both the API and the saved
    packaged data.
5.  It also provides automatic local caching, thanks to
    [BiocFileCache](https://bioconductor.org/packages/BiocFileCache/).

## MultiAssayExperiment Cheatsheet

<a href="https://github.com/waldronlab/cheatsheets/blob/master/MultiAssayExperiment_QuickRef.pdf">
<img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/master/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/>
</a>

## Quick Start

### Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("cBioPortalData", quietly = TRUE))
    BiocManager::install("waldronlab/cBioPortalData")

library(cBioPortalData)
```

## Note

`cBioPortalData` is a work in progress due to changes in data curation
and cBioPortal API specification. Users can view the
`data(studiesTable)` dataset to get an overview of the studies that are
available and currently building as `MultiAssayExperiment`
representations. About 98 % of the studies via the API (`api_build`) and
73 % of the package studies (`pack_build`) are building, these include
additional datasets that were not previously available. Feel free to
file an issue to request prioritization of fixing any of the remaining
datasets.

``` r
cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
```

``` r
table(studies$api_build)
#> 
#> FALSE  TRUE 
#>     5   308

table(studies$pack_build)
#> 
#> FALSE  TRUE 
#>    86   227
```

### API Service

Flexible and granular access to cBioPortal data from
`cbioportal.org/api`. This option is best used with a particular gene
panel of interest. It allows users to download sections of the data with
molecular profile and gene panel combinations within a study.

``` r
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
#>  [2] gbm_tcga_mrna: SummarizedExperiment with 336 rows and 401 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

### Packaged Data Service

This function will download a dataset from the `cbioportal.org/datasets`
website as a packaged tarball and serve it to users as a
`MultiAssayExperiment` object. This option is good for users who are
interested in obtaining all the data for a particular study.

``` r
acc <- cBioDataPack("acc_tcga")
```

``` r
acc
#> A MultiAssayExperiment object of 11 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 11:
#>  [1] cna_hg19.seg: RaggedExperiment with 16080 rows and 90 columns
#>  [2] CNA: SummarizedExperiment with 24776 rows and 90 columns
#>  [3] linear_CNA: SummarizedExperiment with 24776 rows and 90 columns
#>  [4] methylation_hm450: SummarizedExperiment with 15755 rows and 80 columns
#>  [5] mutations_extended: RaggedExperiment with 20166 rows and 90 columns
#>  [6] mutations_mskcc: RaggedExperiment with 20166 rows and 90 columns
#>  [7] RNA_Seq_v2_expression_median: SummarizedExperiment with 20531 rows and 79 columns
#>  [8] RNA_Seq_v2_mRNA_median_all_sample_Zscores: SummarizedExperiment with 20531 rows and 79 columns
#>  [9] RNA_Seq_v2_mRNA_median_Zscores: SummarizedExperiment with 20440 rows and 79 columns
#>  [10] rppa_Zscores: SummarizedExperiment with 191 rows and 46 columns
#>  [11] rppa: SummarizedExperiment with 192 rows and 46 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
