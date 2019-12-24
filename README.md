
# cBioPortalData

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
4.  It provides automatic local caching, thanks to
BiocFileCache.

## MultiAssayExperiment Cheatsheet

<a href="https://github.com/waldronlab/cheatsheets/blob/master/MultiAssayExperiment_QuickRef.pdf"><img src="https://raw.githubusercontent.com/waldronlab/cheatsheets/master/pngs/MultiAssayExperiment_QuickRef.png" width="989" height="1091"/></a>

## Minor note

It is a work in progress, and due to some variation in their formats,
does not yet work for all 268 (as of Dec 2019) datasets. At the time of
writing, it successfully imports 93% of 200 randomly sampled datasets.
Please feel free to file an issue to request prioritization of fixing
any of the remaining datasets.

## Quick Start

### Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("cBioPortalData", quietly = TRUE))
    BiocManager::install("waldronlab/cBioPortalData")

library(cBioPortalData)
```

### API Service

Flexible and granular access to cBioPortal data from
`cbioportal.org/api`.

``` r
cbio <- cBioPortal()
gbm <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "gbm_tcga",
    genePanelId = "IMPACT341",
    molecularProfileIds = c("gbm_tcga_rppa", "gbm_tcga_mrna")
)
#> harmonizing input:
#>   removing 111 colData rownames not in sampleMap 'primary'
gbm
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] gbm_tcga_rppa: SummarizedExperiment with 67 rows and 244 columns
#>  [2] gbm_tcga_mrna: SummarizedExperiment with 334 rows and 401 columns
#> Features:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DFrame
#>  sampleMap() - the sample availability DFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
```

### Packaged data service (legacy)

This function will download a dataset from the `cbioportal.org/datasets`
website as a packaged tarball and serve it to users as a
`MultiAssayExperiment` object.

``` r
laml <- cBioDataPack("laml_tcga")
#> Study file in cache: laml_tcga
#> Working on: /tmp/RtmpXe2G8L/data_CNA.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_RNA_Seq_expression_median.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_RNA_Seq_mRNA_median_Zscores.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_RNA_Seq_v2_expression_median.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_RNA_Seq_v2_mRNA_median_Zscores.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_cna_hg19.seg
#> Parsed with column specification:
#> cols(
#>   ID = col_character(),
#>   chrom = col_double(),
#>   loc.start = col_double(),
#>   loc.end = col_double(),
#>   num.mark = col_double(),
#>   seg.mean = col_double()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_linear_CNA.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_methylation_hm27.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_methylation_hm450.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_double(),
#>   Hugo_Symbol = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   Hugo_Symbol = col_character()
#> )
#> Working on: /tmp/RtmpXe2G8L/data_mutations_extended.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   Entrez_Gene_Id = col_double(),
#>   Start_Position = col_double(),
#>   End_Position = col_double(),
#>   dbSNP_Val_Status = col_logical(),
#>   Score = col_double(),
#>   t_ref_count = col_double(),
#>   t_alt_count = col_double(),
#>   n_ref_count = col_double(),
#>   n_alt_count = col_double(),
#>   Protein_position = col_double(),
#>   Hotspot = col_double(),
#>   RNAVAF_WU = col_double(),
#>   RNAVarReads_WU = col_double(),
#>   stop = col_double(),
#>   NormalVAF_WU = col_double(),
#>   start = col_double(),
#>   TumorVAF_WU = col_double(),
#>   RNARefReads_WU = col_double()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   .default = col_character()
#> )
#> See spec(...) for full column specifications.
#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously

#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously

#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously
#> Working on: /tmp/RtmpXe2G8L/data_mutations_mskcc.txt
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   Entrez_Gene_Id = col_double(),
#>   Start_Position = col_double(),
#>   End_Position = col_double(),
#>   dbSNP_Val_Status = col_logical(),
#>   Score = col_double(),
#>   t_ref_count = col_double(),
#>   t_alt_count = col_double(),
#>   n_ref_count = col_double(),
#>   n_alt_count = col_double(),
#>   Protein_position = col_double(),
#>   Hotspot = col_double(),
#>   RNAVAF_WU = col_double(),
#>   RNAVarReads_WU = col_double(),
#>   stop = col_double(),
#>   NormalVAF_WU = col_double(),
#>   start = col_double(),
#>   TumorVAF_WU = col_double(),
#>   RNARefReads_WU = col_double()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   .default = col_character()
#> )
#> See spec(...) for full column specifications.
#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously

#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously

#> Warning in .find_seqnames_col(df_colnames0, seqnames.field0, xfix): cannnot
#> determine seqnames column unambiguously
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   INITIAL_PATHOLOGIC_DX_YEAR = col_double(),
#>   AGE = col_double(),
#>   PLATELET_COUNT_PRERESECTION = col_double(),
#>   BLAST_COUNT = col_double(),
#>   BASOPHILS_COUNT = col_double(),
#>   ABNORMAL_LYMPHOCYTE_PERCENT = col_double(),
#>   OS_MONTHS = col_double()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   SAMPLE_TYPE_ID = col_double()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   SAMPLE_TYPE_ID = col_double()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   .default = col_character()
#> )
#> See spec(...) for full column specifications.
#> Parsed with column specification:
#> cols(
#>   index = col_double(),
#>   chromosome = col_double(),
#>   region_start = col_double(),
#>   region_end = col_double(),
#>   peak_start = col_double(),
#>   peak_end = col_double(),
#>   enlarged_peak_start = col_double(),
#>   enlarged_peak_end = col_double(),
#>   n_genes_in_region = col_double(),
#>   genes_in_region = col_character(),
#>   n_genes_in_peak = col_double(),
#>   genes_in_peak = col_character(),
#>   n_genes_on_chip = col_character(),
#>   genes_on_chip = col_character(),
#>   `top 3` = col_character(),
#>   amp = col_double(),
#>   cytoband = col_character(),
#>   q_value = col_double()
#> )
#> Warning in .find_with_xfix(df_colnames, prefixes1, prefixes2, start.field, :
#> Multiple prefixes found, using keyword 'region' or taking first one
#> Parsed with column specification:
#> cols(
#>   index = col_double(),
#>   chromosome = col_double(),
#>   region_start = col_double(),
#>   region_end = col_double(),
#>   peak_start = col_double(),
#>   peak_end = col_double(),
#>   enlarged_peak_start = col_double(),
#>   enlarged_peak_end = col_double(),
#>   n_genes_in_region = col_double(),
#>   genes_in_region = col_character(),
#>   n_genes_in_peak = col_double(),
#>   genes_in_peak = col_character(),
#>   n_genes_on_chip = col_character(),
#>   genes_on_chip = col_character(),
#>   `top 3` = col_character(),
#>   amp = col_double(),
#>   cytoband = col_character(),
#>   q_value = col_double()
#> )
#> Warning in .find_with_xfix(df_colnames, prefixes1, prefixes2, start.field, :
#> Multiple prefixes found, using keyword 'region' or taking first one
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
#> Features:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DFrame
#>  sampleMap() - the sample availability DFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
```
