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
        config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L,
            http_version = 0L),
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
    mols <- cbio$getAllMolecularProfilesInStudyUsingGET(studyId = studyId)
    cmols <- httr::content(mols)
    dplyr::bind_rows(cmols)
}

#' Produce small dataset of molecular profile data
#'
#' This function will query the `fetchAllMolecularDataInMolecularProfileUsingPOST`
#' endpoint to obtain data
#'
#' @inheritParams getStudies
#' @param profileId A single string indicating molecular profile ID
#' @param entrezGeneIds A numeric vector indicating entrez gene IDs
#' @param sampleIds A character vector for TCGA sample identifiers
#'
#' @examples
#'
#' cc <- cBioPortal()
#' molecularSlice(cc)
#'
#' @export
molecularSlice <- function(cbio, profileId = "acc_tcga_rna_seq_v2_mrna",
    entrezGeneIds = c(1, 2),
    sampleIds = c("TCGA-OR-A5J1-01",  "TCGA-OR-A5J2-01"))
{
    mols <- cbio$fetchAllMolecularDataInMolecularProfileUsingPOST(
        molecularProfileId = profileId,
        entrezGeneIds = entrezGeneIds,
        sampleIds = sampleIds
    )
    cmols <- httr::content(mols)
    byGene <- dplyr::bind_rows(cmols)
    tidyr::spread(byGene[, c("entrezGeneId", "sampleId", "value")],
        sampleId, value)
}

#' Get a table of all genes
#'
#' Query the API to get a list of Entrez Gene IDs, Hugo symbols, etc.
#'
#' @inheritParams clinicalData
#'
#' @export
geneTable <- function(cbio) {
    gres <- cbio$getAllGenesUsingGET()
    glist <- httr::content(gres)
    glix <- lapply(glist, function(x) {
        if (is.null(x[["cytoband"]]))
            x[["cytoband"]] <- NA_character_
        if (is.null(x[["length"]]))
            x[["length"]] <- NA_integer_
        if (is.null(x[["chromosome"]]))
            x[["chromosome"]] <- NA_character_
        x
    })
    dplyr::bind_rows(glix)
}

#' Return all samples within sample lists
#'
#' Provide a sampleListId and the function will return a CharacterList of
#' associated TCGA barcodes
#'
#' @param cbio An object of class `cBioPortal`
#' @param sampleListIds A character vector of sampleListId as obtained from
#' sampleLists
#'
#' @export
samplesInSampleLists <- function(cbio, sampleListIds = c("acc_tcga_all")) {
    sampleListIds <- setNames(sampleListIds, sampleListIds)
    cnames <- lapply(sampleListIds, function(x) {
        res <- cbio$getAllSampleIdsInSampleListUsingGET(sampleListId = x)
        unlist(httr::content(res))
    })
    IRanges::CharacterList(cnames)
}

#' Provide a sample lists within a study
#'
#' Given a particular `studyId``, this function will return all the available
#' `sampleListId` identifiers
#'
#' @inheritParams clinicalData
#'
#' @export
sampleLists <- function(cbio, studyId = "acc_tcga") {
    slist <- cbio$getAllSampleListsInStudyUsingGET(studyId = studyId)
    slist <- httr::content(slist)
    vapply(slist, `[[`, character(1L), "sampleListId")
}
