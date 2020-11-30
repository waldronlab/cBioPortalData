# Changes in version 2.4.0

## New features

* `cBioPortalData` now allows for gene inputs as either Entrez IDs or Hugo
symbols (#24, @jucor)
* When `gene` inputs are provided, the `by` argument has to agree with the type
of genes provided (either be `entrezGeneId` or `hugoGeneSymbol`).

# Changes in version 2.2.0

## New features

* `studiesTable` includes additional columns `pack_build` and `api_build` to
indicate to the user which datasets have been successfully built as
`MultiAssayExperiment` objects. Users will be notified when a dataset, reported
as not building, is requested from the `cBioDataPack` function.
* Add `sampleIds` argument to `getDataByGenePanel` as part of cache re-work
* Allow more flexibility in the hostname when accessing the API with
`cBioPortal` (@inodb, #16)
* `cBioDataPack` downloads from a more robust repository (AWS S3; @inodb, #22)
* `removePackCache` and `removeDataCache` now remove data from the user's
cache based on inputs to respective functions (`cBioDataPack` and
`cBioPortalData`)

## Bug fixes and minor improvements

* Attempt to merge additional clinical data files from tarballs in
`cBioDataPack`.
* Switch to using `read.delim` instead of `read_tsv` internally to avoid
assigning `NA` to chromosome column
* Use 'PATIENT_ID' when available to determine if experiment data is provided
in the tarball files.
* Add tests using `testthat`
* Update and include percentages of studies successfully imported using
`cBioDataPack` and `cBioPortalData` in the documentation
* Fix read-in when identifiers are numeric instead of character (@jucor, #27)
* Include pagination parameters in `geneTable` function (@xinwei-sher, #29)

# Changes in version 2.0.0

## New features

* Bioconductor release!
* Updated the `README.md` file from R Markdown file.
* Uses the latest version of `rapiclient` on CRAN
* Prepare package for Bioconductor submission
* Include protein metadata as a `RaggedExperiment` from mutation molecular
profiles (TCGA only)

## Bug fixes and minor improvements

* API authentication option removed and not needed

# Changes in version 1.0.1

## New features

* Package supports nearly all study identifiers based on recent tests
* Only a handful of study identifiers are unsuccessful (create an issue to
prioritize).

## Bug fixes and minor improvements

* Make better use of the API return values to craft the sample map for
`MultiAssayExperiment` creation
* Additional data included in the metadata slot of the `MultiAssayExperiment`
object. Future revisions will include this data as `rowData`.
* Change vignette titles for build

# Changes in version 0.1.0

## New features

* `cBioDataPack` allows users to download packaged data objects from
download.cbioportal.org/
* Data packs are cached using `BiocFileCache` to avoid re-downloading
* `cBioPortalData` lets users query the cbioportal.org API and retrieve slices
of data according to gene, molecular profile identifiers, etc.
* Queries through `cBioPortalData` use a caching mechanism to avoid repeat
downloads of data and improve load times
* Both functions return a `MultiAssayExperiment` as the primary data
representation
* Only a number of study datasets are currently possible to load. Issues
can arise with mismatched or munged identifiers
* The cBioPortal API representation is handled by the `AnVIL` package
which makes use of `rapiclient` to provide an automatic R interface to the API

## Bug fixes and minor improvements

* Data pack downloads use an alternative method for download when a download
fails
