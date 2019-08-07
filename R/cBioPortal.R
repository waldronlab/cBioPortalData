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

## curl -X POST "http://www.cbioportal.org/api/gene-panel-data/fetch" -H "accept: application/json" -H "Content-Type: application/json" -d "{ \"sampleMolecularIdentifiers\": [ { \"molecularProfileId\": \"acc_tcga_linear_CNA\", \"sampleId\": \"TCGA-OR-A5J1-01\" } ]}"

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
#' @examples
#' samplesInSampleLists(cb, sampleLists(cb))
#'
#' @export
samplesInSampleLists <- function(cbio, sampleListIds = c("acc_tcga_all")) {
    sampleListIds <- setNames(sampleListIds, sampleListIds)
    cnames <- lapply(sampleListIds, function(x) {
        res <- cbio$getSampleListUsingGET(sampleListId = x)

        res2 <- httr::content(res)
        meta <- res2[names(res2) != "sampleIds"]
        ids <- unlist(res2[["sampleIds"]])
        list(ids = ids, metadata = meta)
    })
    res <- IRanges::CharacterList(cnames[["ids"]])
    metadata(res) <- cnames[["metadata"]]
    res
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
    dplyr::bind_rows(slist)
}

#' @export
genePanel <- function(cbio, panelId = "NSCLC_UNITO_2016_PANEL") {
   gp <- cbio$getGenePanelUsingGET(genePanelId = panelId)
   gp <- httr::content(gp)
   dplyr::bind_rows(gp$genes)
}

# genePanelMolecular <- function(cbio,
#     molecularProfileId = "nsclc_unito_2016_mutations") {
#
#     (molprof <- molecularProfiles(cbio)$molecularProfileId)
#     (samplist <- sampleLists(cbio)$sampleListId)
# #    molprof <- "acc_tcga_linear_CNA"
#     cbio$getGenePanelDataUsingPOST(
#         list(
#             sampleMolecularIdentifiers = data.frame(
#                 molecularProfileId = "acc_tcga_rppa",
#                 sampleListId = "acc_tcga_rppa"
#             )
#         )
#     )
#     getGP <- cbio$getGenePanelDataUsingPOST
#
#     samps <- head(unlist(samplesInSampleLists(cbio, molprof), use.names = FALSE))
#
#     attributes(getGP) <- NULL
#     debugonce(getGP)
#     getGP("acc_tcga_linear_CNA", sampleListId = "acc_tcga_cnaseq")
#
#     ## works !
#     cbio$fetchGenePanelDataInMultipleMolecularProfilesUsingPOST
#     httr:::content(
#         cbio$fetchGenePanelDataInMultipleMolecularProfilesUsingPOST(
#             sampleMolecularIdentifiers =
#                 list( sampleMolecularIdentifiers =
#                     data.frame(
#                         molecularProfileId = "acc_tcga_linear_CNA",
#                         sampleId = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#                     )
#                 )
#         )
#     )
# }
