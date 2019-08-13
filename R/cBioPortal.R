#' @export
cbioportal <- NULL

#' @export
.cBioPortal <- setClass("cBioPortal", contains = "Service")

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
    .cBioPortal(
        Service(
            service = "cBioPortal",
            host = "www.cbioportal.org",
            config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L,
                http_version = 0L),
            package = "cBioPortalData",
            schemes = "http"
        )
    )
}

.invoke_fun <- function(api, name, ...) {
    ops <- names(operations(api))
    if (!name %in% ops)
        stop("<internal> operation name not found in API")
    do.call(`$`, list(api, name))(...)
}


.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_bind <- function(api, name, ...) {
    .bind_content(.invoke_fun(api, name, ...))
}

#' Obtain a table of studies and associated metadata
#'
#' @param cbio An object of class `cBioPortal`
#'
#' @export
getStudies <- function(cbio) {
    query <- .invoke_fun(cbio, "getAllStudiesUsingGET")
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
    dfclin <- .invoke_bind(cbio, "fetchAllClinicalDataInStudyUsingPOST",
        studyId = studyId)
    tidyr::spread(dfclin, clinicalAttributeId, value)
}

#' Produce molecular profiles dataset
#'
#' @inheritParams clinicalData
#'
#' @export
molecularProfiles <- function(cbio, studyId = "acc_tcga",
    projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    projection <- match.arg(projection)
    mols <- .invoke_fun(cbio, "getAllMolecularProfilesInStudyUsingGET",
        studyId = studyId, projection = projection)
    cmols <- httr::content(mols)
    if (projection %in% c("SUMMARY", "ID"))
        dplyr::bind_rows(cmols)
    else
        cmols
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
    byGene <- .invoke_bind(cbio,
        "fetchAllMolecularDataInMolecularProfileUsingPOST",
        molecularProfileId = profileId,
        entrezGeneIds = entrezGeneIds,
        sampleIds = sampleIds
    )
    tidyr::spread(byGene[, c("entrezGeneId", "sampleId", "value")],
        sampleId, value)
}

#' Search through API operations
#'
#' Use a string to filter through all available operations
#'
#' @param cbio An object of class `cBioPortal`
#' @param keyword A string for searching through available operations
#'
#' @examples
#'
#' searchOps(cbio, "mol")
#'
#' @export
searchOps <- function(cbio, keyword) {
    grep(keyword, names(operations(cbio)), value = TRUE, ignore.case = TRUE)
}

#' Get a table of all genes
#'
#' Query the API to get a list of Entrez Gene IDs, Hugo symbols, etc.
#'
#' @inheritParams clinicalData
#'
#' @export
geneTable <- function(cbio, ...) {
    gres <- .invoke_fun(cbio, "getAllGenesUsingGET", ...)
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
#' @examples
#' samplesInSampleLists(cc, "acc_tcga_rppa")
#'
#' @export
samplesInSampleLists <- function(cbio, sampleListIds = c("acc_tcga_all")) {
    sampleListIds <- setNames(sampleListIds, sampleListIds)
    meta <- structure(vector("list", length(sampleListIds)), .Names = sampleListIds)
    res <- lapply(sampleListIds, function(x) {
        res <- cbio$getSampleListUsingGET(sampleListId = x)
        res2 <- httr::content(res)
        meta[[x]] <<- res2[names(res2) != "sampleIds"]
        unlist(res2[["sampleIds"]])
    })
    res <- IRanges::CharacterList(res)
    meta <- dplyr::bind_rows(meta)
    metadata(res) <- meta
    res
}

#' Provide sample lists for a particular study
#'
#' Given a particular `studyId``, this function will return all the available
#' `sampleListId` identifiers
#'
#' @inheritParams clinicalData
#'
#' @examples
#' sampleLists(cc, "acc_tcga")
#'
#' @export
sampleLists <- function(cbio, studyId = "acc_tcga") {
    .invoke_bind(cbio, "getAllSampleListsInStudyUsingGET", studyId = studyId)
}

#' @export
allSamples <- function(cbio, studyId = "acc_tcga") {
    .invoke_bind(cbio, "getAllSamplesInStudyUsingGET", studyId = studyId)
}

#' @export
genePanels <- function(cbio) {
    .invoke_bind(cbio, "getAllGenePanelsUsingGET", studyId = studyId)
}

#' @export
getGenePanel <- function(cbio, panelId = "NSCLC_UNITO_2016_PANEL") {
    .invoke_bind(cbio, "getGenePanelUsingGET", genePanelId = panelId)
}

#' @export
genePanelMolecular <-
    function(cbio, molecularProfileId = "acc_tcga_linear_CNA",
        sampleListId = "acc_tcga_cna", sampleIds = NULL)
{
    .invoke_bind(cbio, "getGenePanelDataUsingPOST",
        molecularProfileId = molecularProfileId,
        sampleListId = list(sampleListId = sampleListId)
    )
}

#' @export
getGenePanelMolecular <-
    function(cbio, molecularProfileId = "acc_tcga_linear_CNA",
        sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01"))
{
    .invoke_bind(cbio,
        "fetchGenePanelDataInMultipleMolecularProfilesUsingPOST",
        sampleMolecularIdentifiers = list(sampleMolecularIdentifiers =
            data.frame(
                molecularProfileId = molecularProfileId,
                sampleId = sampleIds
            )
        )
    )
}

#' @export
getSampleInfo <-
    function(cbio, studyId = "acc_tcga", sampleListIds = NULL,
        projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    projection <- match.arg(projection)
    if (!is.null(sampleListIds))
        queryobj <- list(sampleListIds = sampleListIds)
    else
        queryobj <- list(sampleIdentifiers =
            as.data.frame(
                allSamples(cbio, studyId)[, c("sampleId", "studyId")]
            )
        )

    .invoke_bind(cbio, "fetchSamplesUsingPOST",
        projection = projection, sampleIdentifiers = queryobj
    )
}

