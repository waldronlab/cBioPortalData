.get_cache <- function() {
    cache <- getOption("cBioCache", setCache(verbose = FALSE))

    BiocFileCache(cache)
}

.cache_exists <- function(bfc, rname) {
    file.exists(bfcrpath(bfc, rname, exact = TRUE))
}

.checkSize <- function(cancer_study_id) {

    bfc <- .get_cache()
    study_file <- bfcquery(
        bfc, cancer_study_id, "rname", exact = TRUE)$rpath

    URL <- paste0("http://download.cbioportal.org/", cancer_study_id, ".tar.gz")

    header <- httr::HEAD(URL)$headers
    header_bytes <- as.numeric(header$`content-length`)

    local_bytes <- file.size(study_file)

    message("url: ", header_bytes, " vs. local: ", local_bytes)

    identical(header_bytes, local_bytes)
}

#' @name cBioCache
#'
#' @title Manage cache / download directories for study data
#'
#' @description Managing data downloads is important to save disk space and
#' re-downloading data files. This can be done effortlessly via the integrated
#' BiocFileCache system.
#'
#' @section cBioCache:
#' Get the directory location of the cache. It will prompt the user to create
#' a cache if not already created. A specific directory can be used via
#' `setCache`.
#'
#' @section setCache:
#' Specify the directory location of the data cache. By default, it will
#' go to the user directory as given by:
#' \preformatted{
#'     tools::R_user_dir("cBioPortalData", "cache")
#' }
#'
#' @section removePackCache:
#' Some files may become corrupt when downloading, this function allows
#' the user to delete the tarball associated with a `cancer_study_id` in the
#' cache. This only works for the `cBioDataPack` function. To remove the entire
#' `cBioPortalData` cache, run `unlink("~/.cache/cBioPortalData")`.
#'
#' @param directory The file location where the cache is located. Once set
#' future downloads will go to this folder.
#'
#' @param verbose Whether to print descriptive messages
#'
#' @param ask logical (default TRUE when interactive session) Confirm the file
#' location of the cache directory
#'
#' @param cancer_study_id A single string from `studiesTable` associated
#' with a study tarball
#'
#' @param dry.run logical Whether or not to remove cache files (default TRUE).
#'
#' @param ... For `cBioCache`, arguments passed to `setCache`
#'
#' @md
#'
#' @examples
#'
#' cBioCache()
#'
#' removePackCache("acc_tcga", dry.run = TRUE)
#'
#' @return cBioCache: The path to the cache location
#' @export
cBioCache <- function(...) {
    getOption("cBioCache", setCache(..., verbose = FALSE))
}

#' @rdname cBioCache
#' @export
setCache <-
    function(directory = tools::R_user_dir("cBioPortalData", "cache"),
        verbose = TRUE,
        ask = interactive())
{
    stopifnot(is.character(directory),
        isSingleString(directory), !is.na(directory))

    if (!dir.exists(directory)) {
        if (ask) {
            qtxt <- sprintf(
                "Create cBioPortalData cache at \n    %s? [y/n]: ",
                directory
            )
            answer <- .getAnswer(qtxt, allowed = c("y", "Y", "n", "N"))
            if ("n" == answer)
                stop("'cBioCache' directory not created. Use 'setCache'")
        }
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    options("cBioCache" = directory)

    if (verbose)
        message("cBioPortalData cache directory set to:\n    ",
            directory)
    invisible(directory)
}

#' @rdname cBioCache
#' @export
removePackCache <- function(cancer_study_id, dry.run = TRUE) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (!length(rid)) {
        message("No record found: ", cancer_study_id, ".tar.gz")
    } else if (dry.run) {
            bfcinfo(bfc, rid)
    } else {
        bfcremove(bfc, rid)
        message("Cache record: ", cancer_study_id, ".tar.gz removed")
    }
}

.getHashCache <- function(hashtag) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, hashtag, "rname", exact = TRUE)$rid
    if (!length(rid))
        bfcnew(bfc, hashtag, ext = ".rda")
    else
        bfcquery(bfc, hashtag, "rname", exact = TRUE)$rpath
}

.molDataCache <-
    function(api, studyId = NA_character_, genePanelId = NA_character_,
    molecularProfileIds = NULL, sampleListId = NULL, sampleIds = NULL)
{
    panel <- getGenePanel(api, genePanelId = genePanelId)
    digi <- digest::digest(
        list("getDataByGenePanel", api, studyId, panel,
            sampleIds, molecularProfileIds)
    )
    .getHashCache(digi)
}

.clinDataCache <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    studyId <- force(studyId)
    digi <- digest::digest(list("clinicalData", api, studyId))
    .getHashCache(digi)
}

#' @rdname cBioCache
#'
#' @inheritParams cBioPortalData
#'
#' @inheritParams cBioPortal
#'
#' @examples
#'
#' cbio <- cBioPortal()
#'
#' cBioPortalData(
#'     cbio, by = "hugoGeneSymbol",
#'     studyId = "acc_tcga",
#'     genePanelId = "AmpliSeq",
#'     molecularProfileIds =
#'         c("acc_tcga_rppa", "acc_tcga_linear_CNA", "acc_tcga_mutations")
#' )
#'
#' removeDataCache(
#'     cbio,
#'     studyId = "acc_tcga",
#'     genePanelId = "AmpliSeq",
#'     molecularProfileIds =
#'         c("acc_tcga_rppa", "acc_tcga_linear_CNA", "acc_tcga_mutations"),
#'     dry.run = TRUE
#' )
#'
#' @export
removeDataCache <- function(api, studyId = NA_character_,
    genePanelId = NA_character_, molecularProfileIds = NULL,
    sampleListId = NULL, sampleIds = NULL, dry.run = TRUE, ...)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    formals <- formals()
    call <- std.args(match.call(), formals)
    exargs <- match.args(.portalExperiments, call)
    exargs <- eval.args(exargs)
    exargs <- update.args(exargs)

    cachelocs <- c(experiment_cache = do.call(.molDataCache, exargs),
    clinical_cache = .clinDataCache(exargs[["api"]], exargs[["studyId"]]))

    if (!dry.run)
        vapply(cachelocs, file.remove, logical(1L))

    cachelocs
}
