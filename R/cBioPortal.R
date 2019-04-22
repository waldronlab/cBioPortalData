#' @export
cbioportal <- NULL

#' API Entry function for the cBioPortal data service
#'
#' This function allows the use of the cBioPortal API
#'
#' @return An object of class 'cBioPortal'
#'
#' @importFrom AnVIL Service
#'
#' @export
cBioPortal <- function() {
    Service(
        service = "cBioPortal",
        host = "www.cbioportal.org",
        config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L),
        authenticate_config = FALSE,
        package = "cBioPortalData"
    )
}

#' Obtain a table of studies and associated metadata
#'
#' @param cbio An object of class `cBioPortal`
#'
#' @export
getStudies <- function(cbio) {
    query <- cbio$getAllStudiesUsingGET()
    studies <- httr::content(query)
    studies <- lapply(studies, function(x) {
        if (is.null(x[["pmid"]]))
            x[["pmid"]] <- NA_character_
        if (is.null(x[["citation"]]))
            x[["citation"]] <- NA_character_
        x
    })
    dplyr::bind_rows(studies)
}

#' Obtain clinical data
#'
#' @param cbio An object of class `cBioPortal`
#'
#' @param studyId A single string indicating the "studyId" as taken from
#'     `getStudies`
#'
#' @export
clinicalData <- function(cbio, studyId = "acc_tcga") {
    clinpost <- cbio$fetchAllClinicalDataInStudyUsingPOST(studyId = studyId)
    clindat <- httr::content(clinpost)
    dfclin <- dplyr::bind_rows(clindat)
    tidyr::spread(dfclin, clinicalAttributeId, value)
}

#' Produce molecular profiles dataset
#'
#' @inheritParams clinicalData
#'
#' @export
molecularProfiles <- function(cbio, studyId = "acc_tcga") {
    mols <- cbio$getAllMolecularProfilesInStudyUsingGET(studyId = "acc_tcga")
    cmols <- httr::content(mols)
    dplyr::bind_rows(cmols)
}
