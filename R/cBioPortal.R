## construct a singleton instance for this service

#' @rdname Service
#'
#' @return 'cBioPortal' represents the API of the cBioPortal database
#'
#' @export
cBioPortal <- NULL # assigned in .onLoad, when credentials are available(?)

cBioPortal <-
    function()
{
    Service(
        "cBioPortal",
        host = "https://www.cbioportal.org/api/",
        config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)
    )
}
