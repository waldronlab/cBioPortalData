.get_cache <- function() {
    cache <- getOption("cBioCache", setCache(verbose = FALSE))

    BiocFileCache::BiocFileCache(cache)
}

.cache_exists <- function(bfc, rname) {
    file.exists(bfcrpath(bfc, rname, exact = TRUE))
}

.checkSize <- function(cancer_study_id) {

    bfc <- .get_cache()
    study_file <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rpath

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
#' got to the user's home/.cache and "appname" directory as specified by
#' \link{user_cache_dir}. (default appname: cBioPortalData)
#'
#' @section removeCache:
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
#' @param ... For `cBioCache`, arguments passed to `setCache`
#'
#' @md
#'
#' @examples
#'
#' (cacheloc <- cBioCache())
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
removeCache <- function(cancer_study_id) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (length(rid)) {
        bfcremove(bfc, rid)
        message("Cache record: ", cancer_study_id, ".tar.gz removed")
    } else
        message("No record found: ", cancer_study_id, ".tar.gz")
}

.getHashCache <- function(hashtag) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, hashtag, "rname", exact = TRUE)$rid
    if (!length(rid))
        BiocFileCache::bfcnew(bfc, hashtag, ext = ".rda")
    else
        BiocFileCache::bfcquery(bfc, hashtag, "rname", exact = TRUE)$rpath
}
