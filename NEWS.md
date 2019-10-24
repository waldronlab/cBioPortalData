## CHANGES IN VERSION 0.1.0

### New features

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

### Bug fixes and minor improvements

* Data pack downloads use an alternative method for download when a download
fails
