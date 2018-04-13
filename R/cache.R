.get_cache <- function() {
    cache <- getOption("bio_cache", setCache(verbose = FALSE))

    BiocFileCache::BiocFileCache(cache)
}

.cache_exists <- function(bfc, rname) {
    file.exists(bfcrpath(bfc, rname))
}

#' @name bio_cache
#'
#' @title Manage cache / download directories for study data
#'
#' @description Managing data downloads is important to save disk space and
#' re-downloading data files. This can be done effortlessly via the integrated
#' BiocFileCache system.
#'
#' @section bio_cache:
#' Get the directory location of the cache. It will prompt the user to create
#' a cache if not already created. A specific directory can be used via
#' \code{setCache}.
#'
#' @section setCache:
#' Specify the directory location of the data cache. By default, it will
#' got to the user's home/.cache and "appname" directory as specified by
#' \link{user_cache_dir}. (default appname: MultiAssayExperimentData)
#'
#' @section removeCache:
#' Some files may become corrupt when downloading, this function allows
#' the user to delete the tarball associated with a `cancer_study_id` in the
#' cache.
#'
#' @param directory The file location where the cache is located. Once set
#' future downloads will go to this folder.
#' @param verbose Whether to print descriptive messages
#' @param ask logical (default TRUE when interactive session) Confirm the file
#' location of the cache directory
#' @param cancer_study_id A single string from `studiesTable` associated
#' with a study tarball
#'
#' @md
#'
#' @export
bio_cache <- function(...) {
    getOption("bio_cache", setCache(..., verbose = FALSE))
}

#' @rdname bio_cache
#' @export
setCache <-
function(directory = rappdirs::user_cache_dir("MultiAssayExperimentData"),
    verbose = TRUE,
    ask = interactive())
{
    stopifnot(is.character(directory),
        isSingleString(directory), !is.na(directory))

    if (!dir.exists(directory)) {
        if (ask) {
            qtxt <- sprintf(
                "Create MultiAssayExperimentData cache at \n    %s? [y/n]: ",
                directory
            )
            answer <- .getAnswer(qtxt, allowed = c("y", "Y", "n", "N"))
            if ("n" == answer)
                stop("'bio_cache' directory not created. Use 'setCache'")
        }
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    options("bio_cache" = directory)

    if (verbose)
        message("MultiAssayExperimentData cache directory set to:\n    ",
            directory)
    invisible(directory)
}

#' @rdname bio_cache
#' @export
removeCache <- function(cancer_study_id) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname")$rid
    if (length(rid)) {
        bfcremove(bfc, rid)
        message("Cache record: ", cancer_study_id, ".tar.gz removed")
    } else
        message("No record found: ", cancer_study_id, ".tar.gz")
}
