#' [Defunct] Available studies from cBioPortal
#'
#' Note. This dataset has been replaced by `getStudies()` from the
#' cBioPortal API.
#'
#' @format A data frame with 220 rows and 4 variables:
#' \describe{
#'     \item{cancer_study_id}{
#'         The study code used for input to `cBioDataPack`
#'     }
#'     \item{study_name}{
#'         A descriptive study title containing data center and year
#'     }
#'     \item{description}{
#'         A longer description of the study
#'     }
#'     \item{URL}{
#'         Associated study URLs
#'     }
#' }
#' @docType data
#' @name studiesTable
#'
#' @author Marcel Ramos \email{marcel.ramos@@roswellpark.org}
#'
#' @references \url{http://www.cbioportal.org/datasets},
#'   \url{https://github.com/cBioPortal/cgdsr}
#' @keywords data
NULL
