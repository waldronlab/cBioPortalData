#' @export
cbioportal <- NULL

#' API Entry function for the cBioPortal data service
#'
#' This function allows the use of the cBioPortal API
#'
#' @return An object of class 'cBioPortal'
#'
#' @export
cBioPortal <- function() {
    AnVIL::Service(
        service = "cBioPortal",
        host = "www.cbioportal.org",
        config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L),
        authenticate_config = FALSE,
        package = "cBioPortalData"
    )
}
