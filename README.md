# cBioPortalData: Datasets using the MultiAssayExperiment data class

This project aims to import all cBioPortal datasets as [MultiAssayExperiment](http://bioconductor.org/packages/MultiAssayExperiment/)
objects in Bioconductor. It offers some advantages over the CDGS-R package:

1. The MultiAssayExperiment class explicitly links all assays to the patient clinical/pathological data
2. The MultiAssayExperiment class provides a [flexible API](https://github.com/waldronlab/MultiAssayExperiment/wiki/MultiAssayExperiment-API) ([cheatsheet](http://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment_cheatsheet.pdf)), including harmonized subsetting and reshaping to convenient wide and long formats.
3. It provides complete datasets, not just for subsets of genes
4. It provides automatic local caching, thanks to BiocFileCache.

It is a work in progress, and due to some variation in their formats, does not yet work for 
all 213 (as of Oct 2018) datasets. At the time of writing, it successfully imports 145 datasets. Please feel free to file an issue 
to request prioritization of fixing any of the remaining datasets.

This package will eventually be published in Bioconductor once all datasets are successfully imported.
