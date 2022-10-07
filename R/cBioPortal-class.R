.api_header <- function(x) x@api_header

#' @name cBioPortal-class
#'
#' @title A class for representing the cBioPortal API protocol
#'
#' @description The `cBioPortal` class is a representation of the `cBioPortal`
#'     API protocol that directly inherits from the `Service` class in the
#'     `AnVIL` package. For more information, see the
#'     \link[AnVIL:Service]{AnVIL} package.
#'
#' @details This class takes the static API as provided at
#'     \url{https://www.cbioportal.org/api/api-docs} and creates an R object
#'     with the help from underlying infrastructure (i.e.,
#'     \link[rapiclient:rapiclient-package]{rapiclient} and
#'     \link[AnVIL:Service]{AnVIL}) to give the user a unified representation
#'     of the API specification provided by the cBioPortal group. Users are not
#'     expected to interact with this class other than to use it as input
#'     to the functionality provided by the rest of the package.
#'
#' @importFrom methods new
#'
#' @return A \code{cBioPortal} class instance
#'
#' @seealso  \link{cBioPortal}, \linkS4class{Service}
#'
#' @md
#'
#' @examples
#'
#' cBioPortal()
#'
#' @exportClass cBioPortal
.cBioPortal <- setClass(
    "cBioPortal",
    contains = "Service",
    slots = c(api_header = "character")
)

#' @describeIn cBioPortal-class
#'
#' @importFrom AnVIL operations
#' @importFrom methods callNextMethod
#'
#' @param x A \linkS4class{Service} instance or API representation as
#'     given by the \link{cBioPortal} function.
#'
#' @inheritParams AnVIL::operations
#'
#' @export
setMethod(
    "operations", "cBioPortal",
    function(x, ..., .deprecated = FALSE)
{
    callNextMethod(
        x, .headers = .api_header(x), ..., .deprecated = .deprecated
    )
})
