#' Delete downloaded study tarballs from cache
#'
#' Some files may become corrupt when downloading, this function allows
#' the user to delete the tarball associated with a `cancer_study_id``
#'
#' @param cancer_study_id A single string from `studiesTable` associated
#' with a study tarball
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
