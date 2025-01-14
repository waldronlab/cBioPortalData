.api_header <- function(x) x@api_header

#' @name cBioPortal-class
#'
#' @title A class for representing the cBioPortal API protocol
#'
#' @description The `cBioPortal` class is a representation of the `cBioPortal`
#'     API protocol that directly inherits from the `Service` class in the
#'     `AnVIL` package. For more information, see the
#'     [AnVIL][AnVIL::Service-class] package.
#'
#' @details This class takes the static API as provided at
#'     <https://www.cbioportal.org/api/v2/api-docs> and creates an R object
#'     with the help from underlying infrastructure (i.e.,
#'     [rapiclient][rapiclient::rapiclient-package] and
#'     [AnVIL][AnVIL::Service-class]) to give the user a unified representation
#'     of the API specification provided by the cBioPortal group. Users are not
#'     expected to interact with this class other than to use it as input
#'     to the functionality provided by the rest of the package.
#'
#' @importFrom methods new
#'
#' @return A `cBioPortal` class instance
#'
#' @seealso  [cBioPortal], [AnVIL][AnVIL::Service-class]
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

#' @describeIn cBioPortal-class List all the `operations` available with the
#'   cBioPortal API object, e.g., `api$operation`
#'
#' @importFrom AnVIL operations
#' @importFrom methods callNextMethod
#'
#' @param x A [AnVIL][AnVIL::Service-class] instance or API representation as
#'     given by the [cBioPortal] function.
#'
#' @inheritParams AnVIL::operations
#'
#' @exportMethod operations
setMethod(
    "operations", "cBioPortal",
    function(x, ..., .deprecated = FALSE)
{
    callNextMethod(
        x, .headers = .api_header(x), ..., .deprecated = .deprecated
    )
})
