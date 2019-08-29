.invoke_fun <- function(api, name, ...) {
    if (!is(api, "cBioPortal"))
        stop("Provide a 'cBioPortal' class API object")
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

#' @export
cbioportal <- NULL

#' @export
.cBioPortal <- setClass("cBioPortal", contains = "Service")

#' @rdname cBioPortal
#'
#' @title The interface to the cBioPortal API Service
#'
#' @description This function allows the use of the cBioPortal API
#'
#' @param cbio An object of class `cBioPortal`
#'
#' @param studyId character(1) Indicates the "studyId" as taken from
#'     `getStudies`
#'
#' @param keyword character(1) Keyword or pattern for searching through
#'     available operations
#'
#' @param molecularProfileId character(1) Indicates a molecular profile ID
#'
#' @param molecularProfileIds character() A vector of molecular profile IDs
#'
#' @param entrezGeneIds numeric() A vector indicating entrez gene IDs
#'
#' @param sampleIds character() TCGA sample identifiers
#'
#' @param genePanelId character(1) Identifies the gene panel, as obtained
#'     from the `genePanels` function
#'
#' @return
#'
#'     cBioPortal: An object of class 'cBioPortal'
#'
#'     cBioPortaldata: An object of class 'MultiAssayExperiment'
#'
#' @importFrom AnVIL Service
#'
#' @examples
#'
#' cc <- cBioPortal()
#' molecularSlice(cc)
#'
#' searchOps(cc, "mol")
#'
#' samplesInSampleLists(cc, "acc_tcga_rppa")
#'
#' @export
cBioPortal <- function() {
    .cBioPortal(
        Service(
            service = "cBioPortal",
            host = "www.cbioportal.org",
            config = httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L,
                http_version = 0L),
            api_url = "https://www.cbioportal.org/api/api-docs",
            package = "cBioPortalData",
            schemes = "http"
        )
    )
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * getStudies - Obtain a table of studies and associated metadata
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

#' @name cBioPortal
#'
#' @section clinicalData:
#'      Obtain clinical data for a particular study identifier
#'
#' @export
clinicalData <- function(cbio, studyId = "acc_tcga") {
    pttable <- .invoke_bind(cbio,
        "getAllPatientsInStudyUsingGET", studyId = studyId)
    ptrow <- lapply(pttable[["patientId"]], function(pt) {
        .invoke_bind(cbio,
            "getAllClinicalDataOfPatientInStudyUsingGET",
            studyId = studyId, patientId = pt)
    })
    clin <- dplyr::bind_rows(ptrow)
    tidyr::spread(clin, clinicalAttributeId, value)
}

#' @name cBioPortal
#'
#' @section Molecular Profiles:
#'      * molecularProfiles - Produce a molecular profiles dataset for a given
#'      study identifier
#'
#' @inheritParams cBioPortal
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

#' @name cBioPortal
#'
#' @section Molecular Profiles:
#'     * molecularSlice - Produce a dataset of molecular profile data based on
#'     `molecularProfileId`, `entrezGeneIds`, and `sampleIds`
#'
#' @export
molecularSlice <- function(cbio,
    molecularProfileId = "acc_tcga_rna_seq_v2_mrna",
    entrezGeneIds = c(1, 2),
    sampleIds = c("TCGA-OR-A5J1-01",  "TCGA-OR-A5J2-01"))
{
    byGene <- .invoke_bind(cbio,
        "fetchAllMolecularDataInMolecularProfileUsingPOST",
        molecularProfileId = molecularProfileId,
        entrezGeneIds = entrezGeneIds,
        sampleIds = sampleIds
    )
    if ("message" %in% names(byGene)) {
        warning(byGene[["message"]])
        dplyr::tibble()
    } else
        tidyr::spread(byGene[, c("entrezGeneId", "sampleId", "value")],
            sampleId, value)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * searchOps - Search through API operations with a keyword
#'
#' @export
searchOps <- function(cbio, keyword) {
    grep(keyword, names(operations(cbio)), value = TRUE, ignore.case = TRUE)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * geneTable - Get a table of all genes by 'entrezGeneId' or
#'     'hugoGeneSymbol'
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

#' @name cBioPortal
#'
#' @section Sample Data:
#'     * samplesInSampleLists - get all samples associated with a 'sampleListId'
#'
#' @param sampleListIds character() A vector of 'sampleListId' as obtained from
#'     `sampleLists`
#'
#' @export
samplesInSampleLists <- function(cbio, sampleListIds = "acc_tcga_all") {
    sampleListIds <- setNames(sampleListIds, sampleListIds)
    meta <- structure(vector("list", length(sampleListIds)),
        .Names = sampleListIds)
    res <- lapply(sampleListIds, function(x) {
        res <- .invoke_fun(cbio, "getSampleListUsingGET", sampleListId = x)
        res2 <- httr::content(res)
        meta[[x]] <<- res2[names(res2) != "sampleIds"]
        unlist(res2[["sampleIds"]])
    })
    res <- IRanges::CharacterList(res)
    meta <- dplyr::bind_rows(meta)
    metadata(res) <- meta
    res
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * sampleLists - obtain all `sampleListIds` for a particular `studyId`
#'
#' @examples
#' sampleLists(cc, "acc_tcga")
#'
#' @export
sampleLists <- function(cbio, studyId = "acc_tcga") {
    .invoke_bind(cbio, "getAllSampleListsInStudyUsingGET", studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * allSamples - obtain all samples within a particular `studyId`
#'
#' @export
allSamples <- function(cbio, studyId = "acc_tcga") {
    .invoke_bind(cbio, "getAllSamplesInStudyUsingGET", studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * genePanels - Show all available gene panels
#'
#' @export
genePanels <- function(cbio) {
    .invoke_bind(cbio, "getAllGenePanelsUsingGET")
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getGenePanels - Obtain the gene panel for a particular 'genePanelId'
#'
#' @export
getGenePanel <- function(cbio, genePanelId = "NSCLC_UNITO_2016_PANEL") {
    res <- .invoke_fun(cbio, "getGenePanelUsingGET", genePanelId = genePanelId)
    res <- httr::content(res)[["genes"]]
    dplyr::bind_rows(res)
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * genePanelMolecular - get gene panel data for a paricular
#'     `molecularProfileId` and `sampleListId` combination
#'
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

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getGenePanelMolecular - get gene panel data for a combination of
#'     `molecularProfileId` and `sampleListId` vectors
#'
#' @export
getGenePanelMolecular <-
    function(cbio, molecularProfileIds =
        c("acc_tcga_linear_CNA", "acc_tcga_rppa"),
        sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01"))
{
    SampMolIds <- S4Vectors::expand.grid(
        molecularProfileId = molecularProfileIds,
        sampleId = sampleIds
    )
    .invoke_bind(cbio,
        "fetchGenePanelDataInMultipleMolecularProfilesUsingPOST",
        sampleMolecularIdentifiers =
            list(sampleMolecularIdentifiers = SampMolIds)
    )
}

#' @name cBioPortal
#'
#' @section Sample Data:
#'     * getSampleInfo - Obtain sample metadata for a particular `studyId` or
#'     `sampleListId`
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

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getDataByGenePanel - Download data for a gene panel and
#'     `molecularProfileId` combination, optionally a `sampleListId` can be
#'     provided.
#'
#' @param by character(1) Whether to use 'entrezGeneId' or 'hugoGeneSymbol'
#'     as row metadata
#'
#' @export
getDataByGenePanel <-
    function(cbio, by = c("entrezGeneId", "hugoGeneSymbol"),
        genePanelId = "bait_v5",
        studyId = "acc_tcga",
        molecularProfileId = "acc_tcga_rna_seq_v2_mrna",
        sampleListId = NULL)
{
    by <- match.arg(by)
    if (!is.null(sampleListId))
        samples <- samplesInSampleLists(cbio, sampleListId)[[1L]]
    else
        samples <- allSamples(cbio, studyId)[["sampleId"]]

    panel <- getGenePanel(cbio, genePanelId = genePanelId)
    molecularData <- molecularSlice(cbio,
        molecularProfileId = molecularProfileId,
        entrezGeneIds = panel[["entrezGeneId"]],
        sampleIds = samples)

    if (identical(by , "hugoGeneSymbol"))
        dplyr::bind_cols(
            hugoGeneSymbol = unlist(
                panel[match(molecularData[["entrezGeneId"]],
                    panel[["entrezGeneId"]]), by]),
            molecularData[, names(molecularData) != "entrezGeneId"]
        )
    else
        molecularData
}

.portalExperiments <-
    function(cbio, by = c("entrezGeneId", "hugoGeneSymbol"),
        genePanelId = "bait_v5",
        studyId = "acc_tcga",
        molecularProfileIds = NULL,
        sampleListId = NULL)
{
    by <- match.arg(by)
    if (is.null(molecularProfileIds))
        molecularProfileIds <-
            molecularProfiles(cbio, studyId)[["molecularProfileId"]]

    molecularProfileIds <- setNames(molecularProfileIds, molecularProfileIds)

    expers <- lapply(molecularProfileIds, function(molprof) {
        moldata <- getDataByGenePanel(cbio, by = by,
            genePanelId = genePanelId, studyId = studyId,
            molecularProfileId = molprof, sampleListId = sampleListId)
        moldata <- as.data.frame(moldata)
        rownames(moldata) <- moldata[[by]]
        moldata <- data.matrix(moldata[, names(moldata) != by])
        SummarizedExperiment(moldata)
    })
    as(Filter(length, expers), "List")
}

#' Download data from the cBioPortal API
#'
#' Obtain a `MultiAssayExperiment` object for a particular gene panel,
#' `studyId`, `molecularProfileIds`, and `sampleListIds` combination. Default
#' `molecularProfileIds` and `sampleListIds` are set to NULL for including all
#' data.
#'
#' @inheritParams cBioPortal
#'
#' @examples
#'
#' cb <- cBioPortal()
#' cBioPortalData(cb, by = "hugoGeneSymbol")
#'
#' @export
cBioPortalData <-
    function(cbio, studyId,
        genePanelId,
        molecularProfileIds = NULL,
        sampleListId = NULL,
        by = c("entrezGeneId", "hugoGeneSymbol")
    )
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    if (missing(studyId))
        stop("Provide a valid 'studyId' from 'getStudies()'")
    if (missing(genePanelId))
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    .checkIdValidity(cbio, element = studyId, ename = "studyId")
    .checkIdValidity(cbio, element = genePanelId, ename = "genePanelId")
    if (!is.null(molecularProfileIds))
        .checkIdValidity(cbio, element = molecularProfileIds,
            ename = "genePanelId")
    if (!is.null(sampleListId))
        .checkIdValidity(cbio, element = sampleListId, ename = "sampleListId")

    explist <- .portalExperiments(cbio = cbio, by = by,
        genePanelId = genePanelId, studyId = studyId,
        molecularProfileIds = molecularProfileIds,
        sampleListId = sampleListId)

    clin <- clinicalData(cbio, studyId = studyId)
    clin <- as.data.frame(clin)
    rownames(clin) <- clin[["patientId"]]
    if (all(startsWith(rownames(clin), "TCGA")))
        idConvert <- TCGAutils::TCGAbarcode
    else
        idConvert <- identical
    sampmap <- TCGAutils::generateMap(experiments = explist,
        colData = clin, idConverter = idConvert)

    MultiAssayExperiment(explist, clin, sampmap)
}
