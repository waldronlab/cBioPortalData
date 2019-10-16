.validStudyID <- function(cancer_study_id) {

    if (missing(cancer_study_id))
        stop("Provide a valid 'cancer_study_id' from 'studiesTable'")

    stopifnot(is.character(cancer_study_id),
        !is.na(cancer_study_id), length(cancer_study_id) == 1L)

    cancer_study_id <- tolower(cancer_study_id)
    ## Load dataset to envir
    loc_data <- new.env(parent = emptyenv())
    data("studiesTable", envir = loc_data)
    studiesTable <- loc_data[["studiesTable"]]

    ## Ensure study ID is valid
    inTable <- cancer_study_id %in% studiesTable[["cancer_study_id"]]

    if (!inTable)
        stop("Study identifier not found in look up table")
    else
        inTable
}

.download_data_file <-
    function(fileURL, cancer_study_id, verbose = FALSE, force = FALSE)
{
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (!length(rid)) {
        rid <- names(bfcadd(bfc, cancer_study_id, fileURL, download = FALSE))
    }
    if (!.cache_exists(bfc, cancer_study_id) || force) {
        if (verbose)
            message("Downloading study file: ", cancer_study_id, ".tar.gz")
            bfcdownload(bfc, rid, ask = FALSE)
    } else
        message("Study file in cache: ", cancer_study_id)

    bfcrpath(bfc, rids = rid)
}

.manageLocalFile <- function(cancer_study_id, inpath) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (!length(rid))
        stop("Can't update non-existing cache item")

    cachedir <- bfccache(bfc)
    finalname <- paste0(gsub("file", "", basename(tempfile())), "_",
        cancer_study_id, ".tar.gz")
    fileLoc <- file.path(cachedir, finalname)
    file.copy(inpath, fileLoc)

    bfcupdate(bfc, rids = rid, rpath = fileLoc)

    file.remove(inpath)

    bfcrpath(bfc, rids = rid)
}

.altDownload <- function(fileURL, cancer_study_id, verbose = FALSE) {
    if (!requireNamespace("curl"))
        stop("Download the 'curl' package to recover from download errors")

    if (verbose)
        message("Downloading study file: ", cancer_study_id, ".tar.gz")

    tmpFile <- file.path(tempdir(), paste0(cancer_study_id, ".tar.gz"))
    download.file(fileURL, destfile = tmpFile, quiet = TRUE, method = "wget")

    .manageLocalFile(cancer_study_id, tmpFile)
}

#' Download and cache study dataset
#'
#' Provide a `cancer_study_id` from the `studiesTable` and retrieve
#' the study tarball from cBioPortal
#'
#' @param cancer_study_id The cBioPortal study identifier as in
#'     \url{https://cbioportal.org/webAPI}
#'
#' @param use_cache logical (default TRUE) create the default cache location
#' and use it to track downloaded data. If data found in the cache, data will
#' not be re-downloaded. A path can also be provided to data cache location.
#'
#' @param force logical (default FALSE) whether to force re-download data from
#' remote location
#'
#' @examples
#'
#' downloadStudy("laml_tcga_pub", use_cache = tempdir())
#'
#' @keywords internal
downloadStudy <- function(cancer_study_id, use_cache = TRUE, force = FALSE)
{
    .validStudyID(cancer_study_id)

    url_location <- "http://download.cbioportal.org"
    url_file <- file.path(url_location, paste0(cancer_study_id, ".tar.gz"))

    if (is.character(use_cache) && length(use_cache) == 1L)
        cbioCache(directory = use_cache)
    else if (isTRUE(use_cache))
        cbioCache()
    else
        stop("Use 'setCache' or specify a download location")

    tryCatch({
        .download_data_file(url_file, cancer_study_id, verbose = TRUE,
            force = force)
        },
        error = function(cond) {
            message("\n", cond)
            message("\nRetrying download with alternative function...")
            .altDownload(url_file, cancer_study_id, verbose = TRUE)
        }
    )
}

