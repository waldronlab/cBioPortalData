## construct a singleton instance for this service

#' @rdname Service
#'
#' @return 'cBioPortal' represents the API of the cBioPortal database
#'
#' @export
cbioportal <- NULL # assigned in .onLoad

## @export
# .cBioPortal <- setClass("cBioPortal", contains = "Service")

#' API Entry function for the cBioPortal data service
#'
#' This function allows the use of the cBioPortal API
#'
#' @return An object of class 'cBioPortal'
#'
#' @export
cBioPortal <-
    function()
{
    # .cBioPortal(
        AnVIL:::Service(
            "cBioPortal",
            host = "www.cbioportal.org",
            config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L),
            authenticate_config = FALSE,
            package = "cBioPortalData"
         )
    # )
}
