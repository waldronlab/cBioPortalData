#' @name cBio_cache
#'
#' @title Manage cache / download directories for study data
#'
#' @description Managing data downloads is important to save disk space and
#' re-downloading data files. This can be done effortlessly via the integrated
#' BiocFileCache system.
#'
#' @section setCache:
#' Specify the directory location of the data cache. By default, it will
#' got to the user's home directory and project directory as specified by
#' \link{user_cache_dir}
#'
#' @section removeCache:
#' Some files may become corrupt when downloading, this function allows
#' the user to delete the tarball associated with a `cancer_study_id``
#'
#' @param cancer_study_id A single string from `studiesTable` associated
#' with a study tarball
#' @param directory The file location where the cache is located. Once set
#' future downloads will go to this folder.
#'
#' @md
#'
#' @export
removeCache <- function(cancer_study_id) {
    bfc <- .get_cache(TRUE)
    rid <- bfcquery(bfc, cancer_study_id, "rname")$rid
    if (length(rid)) {
        bfcremove(bfc, rid)
        message("Cache record: ", cancer_study_id, ".tar.gz removed")
    } else
        message("No record found: ", cancer_study_id, ".tar.gz")
}

#' @rdname cBio_cache
#' @export
setCache <-
function(directory = rappdirs::user_cache_dir("MultiAssayExperimentData"),
    verbose = TRUE,
    ask = interactive())
{
    create_path <- function(directory) {
        dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }

    stopifnot(is.character(directory),
        isSingleString(directory), !is.na(directory))

    if (!ask)
        create_path(directory)
    else {
        qtxt <- sprintf(
            "Create MultiAssayExperimentData cache at \n    %s? [y/n]: ",
            directory
        )
        answer <- .getAnswer(qtxt, allowed = c("y", "Y", "n", "N"))
        if ("n" == answer)
            stop("'cBio_cache' directory will not be created")
    }
    create_path(directory)
    options("cBio_cache" = directory)

    if (verbose)
        message("MultiAssayExperimentData cache directory set to: ",
                directory)
    invisible(directory)
}

#' @rdname cBio_cache
#' @export
cBio_cache <- function() {
    cache_dir <- getOption("cBio_cache", setCache(verbose = FALSE))
    if (!dir.exists(cache_dir))
        setCache(cache_dir)

    return(cache_dir)
}
