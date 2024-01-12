
# cBioPortalData

<!-- start badges here -->

[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/cBioPortalData.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/cBioPortalData)
[![Platforms](http://www.bioconductor.org/shields/availability/release/cBioPortalData.svg)](https://www.bioconductor.org/packages/release/bioc/html/cBioPortalData.html#archives)
[![Downloads](https://www.bioconductor.org/shields/downloads/release/cBioPortalData.svg)](https://bioconductor.org/packages/stats/bioc/cBioPortalData/)
<!-- end badges here -->

## cBioPortal data and MultiAssayExperiment

### Overview

The `cBioPortalData` R package aims to import cBioPortal datasets as
*[MultiAssayExperiment](https://bioconductor.org/packages/3.19/MultiAssayExperiment)*
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
    *[BiocFileCache](https://bioconductor.org/packages/3.19/BiocFileCache)*

## MultiAssayExperiment Cheatsheet

<a href="https://github.com/waldronlab/cheatsheets/blob/master/MultiAssayExperiment_QuickRef.pdf">
<img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/master/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/>
</a>

## Quick Start

### Installation

To install from Bioconductor (recommended for most users, this will
install the release or development version corresponding to your version
of Bioconductor):

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cBioPortalData")
```

Developers may want to install from GitHub for bleeding-edge updates
(although this is generally not necessary because changes here are also
pushed to
[bioc-devel](https://contributions.bioconductor.org/use-devel.html)).
Note that developers must be working with the development version of
Bioconductor; see
[bioc-devel](https://contributions.bioconductor.org/use-devel.html) for
details.

``` r
if (!require("cBioPortalData", quietly = TRUE))
    BiocManager::install("waldronlab/cBioPortalData")
```

To load the package:

``` r
library(cBioPortalData)
```

## Note

`cBioPortalData` is a work in progress due to changes in data curation
and cBioPortal API specification. Users can view the
`data(studiesTable)` dataset to get an overview of the studies that are
available and currently building as `MultiAssayExperiment`
representations. About 89 % of the studies via the API (`api_build`) and
93 % of the package studies (`pack_build`) are building, these include
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
#>    40   340

table(studies$pack_build)
#> 
#> FALSE  TRUE 
#>    28   352
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
#>  [1] gbm_tcga_mrna: SummarizedExperiment with 336 rows and 401 columns
#>  [2] gbm_tcga_rppa: SummarizedExperiment with 67 rows and 244 columns
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
#> A MultiAssayExperiment object of 10 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 10:
#>  [1] cna_hg19.seg: RaggedExperiment with 16080 rows and 90 columns
#>  [2] cna: SummarizedExperiment with 24776 rows and 90 columns
#>  [3] linear_cna: SummarizedExperiment with 24776 rows and 90 columns
#>  [4] methylation_hm450: SummarizedExperiment with 15754 rows and 80 columns
#>  [5] mrna_seq_v2_rsem_zscores_ref_all_samples: SummarizedExperiment with 20531 rows and 79 columns
#>  [6] mrna_seq_v2_rsem_zscores_ref_diploid_samples: SummarizedExperiment with 20440 rows and 79 columns
#>  [7] mrna_seq_v2_rsem: SummarizedExperiment with 20531 rows and 79 columns
#>  [8] mutations: RaggedExperiment with 20166 rows and 90 columns
#>  [9] rppa_zscores: SummarizedExperiment with 191 rows and 46 columns
#>  [10] rppa: SummarizedExperiment with 192 rows and 46 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
