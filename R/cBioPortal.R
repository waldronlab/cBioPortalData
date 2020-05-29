utils::globalVariables(c("clinicalAttributeId", "value", "sampleId"))

#' @name cBioPortal-class
#'
#' @title A class for representing the cBioPortal API protocol
#'
#' @description The `cBioPortal` class is a representation of the cBioPortal
#'     API protocol that directly inherits from the `Service` class in the
#'     `AnVIL` package. For more information, see the 'AnVIL' package.
#'
#' @details This class takes the static API as provided at
#'     \url{https://www.cbioportal.org/api/api-docs} and creates an R object
#'     with the help from underlying infrastructure (i.e., 'rapiclient' and
#'     'AnVIL') to give the user a unified representation of the API
#'     specification provided by the cBioPortal group. Users are not
#'     expected to interact with this class other than to use it as input
#'     to the functionality provided by the rest of the package.
#'
#' @importFrom methods new
#'
#' @seealso  \link{cBioPortal}, \link[AnVIL]{Service}
#'
#' @examples
#'
#' cBioPortal()
#'
#' @export
.cBioPortal <- setClass("cBioPortal", contains = "Service")

#' @rdname cBioPortal
#'
#' @aliases cBioPortal
#'
#' @title The R interface to the cBioPortal API Data Service
#'
#' @description This section of the documentation lists the functions that
#'     allow users to access the cBioPortal API. The main representation of the
#'     API can be obtained from the `cBioPortal` function. The supporting
#'     functions listed here give access to specific parts of the API and
#'     allow the user to explore the API with individual calls. Many of the
#'     functions here are listed for documentation purposes and are
#'     recommended for advanced usage only. Users should only need to use the
#'     `cBioPortalData` main function to obtain data.
#'
#' @param api An API object of class `cBioPortal` from the `cBioPortal`
#'     function
#'
#' @param hostname character(1) The internet location of the service
#'     (default: 'www.cbioportal.org')
#'
#' @param protocol character(1) The internet protocol used to access the
#'     hostname (default: 'https')
#'
#' @param api. character(1) The directory location of the API protocol within
#'     the hostname (default: '/api/api-docs')
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
#' @param sampleIds character() Sample identifiers
#'
#' @param genePanelId character(1) Identifies the gene panel, as obtained
#'     from the `genePanels` function
#'
#' @return
#'
#'     cBioPortal: An API object of class 'cBioPortal'
#'
#'     cBioPortalData: A data object of class 'MultiAssayExperiment'
#'
#' @importFrom AnVIL Service
#'
#' @examples
#' cbio <- cBioPortal()
#'
#' getStudies(api = cbio)
#'
#' searchOps(api = cbio, keyword = "molecular")
#'
#' clinicalData(api = cbio, studyId = "acc_tcga")
#'
#' molecularProfiles(api = cbio, studyId = "acc_tcga")
#'
#' molecularData(
#'     api = cbio,
#'     molecularProfileId = "acc_tcga_rna_seq_v2_mrna",
#'     entrezGeneIds = c(1, 2),
#'     sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#' )
#'
#' sampleLists(api = cbio, studyId = "acc_tcga")
#'
#' samplesInSampleLists(
#'     api = cbio,
#'     sampleListIds = c("acc_tcga_rppa", "acc_tcga_cnaseq")
#' )
#'
#' genePanels(api = cbio)
#'
#' getGenePanel(api = cbio, genePanelId = "IMPACT341")
#'
#' @export
cBioPortal <- function(
    hostname = "www.cbioportal.org", protocol = "https", api. = "/api/api-docs"
) {
    apiUrl <- paste0(protocol, "://", hostname, api.)
    .cBioPortal(
        Service(
            service = "cBioPortal",
            host = hostname,
            config = httr::config(
                ssl_verifypeer = 0L, ssl_verifyhost = 0L, http_version = 0L
            ),
            authenticate = FALSE,
            api_url = apiUrl,
            package = "cBioPortalData",
            schemes = protocol
        )
    )
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * getStudies - Obtain a table of studies and associated metadata
#'
#' @export
getStudies <- function(api) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    digi <- .inputDigest(match.call(), "getStudies")
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
        query <- .invoke_fun(api, "getAllStudiesUsingGET")
        studies <- httr::content(query)
        studies <- lapply(studies, function(x) {
            if (is.null(x[["pmid"]]))
                x[["pmid"]] <- NA_character_
            if (is.null(x[["citation"]]))
                x[["citation"]] <- NA_character_
            x
        })
        save(studies, file = cacheloc)
    }
    dplyr::bind_rows(studies)
}

#' @name cBioPortal
#'
#' @section Patient Data:
#'      * clinicalData - Obtain clinical data for a particular study identifier
#'      ('studyId')
#'
#' @export
clinicalData <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    digi <- .inputDigest(match.call(), "clinicalData")
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
        pttable <- .invoke_bind(
            api = api, name = "getAllPatientsInStudyUsingGET",
            use_cache = TRUE, studyId = studyId
        )
        ptrow <- lapply(pttable[["patientId"]], function(pt) {
            .invoke_bind(
                api = api, name = "getAllClinicalDataOfPatientInStudyUsingGET",
                use_cache = TRUE, studyId = studyId, patientId = pt
            )
        })
        save(ptrow, file = cacheloc)
    }
    clin <- dplyr::bind_rows(ptrow)
    tidyr::pivot_wider(data = clin, names_from = "clinicalAttributeId",
        values_from = "value")
}

#' @name cBioPortal
#'
#' @section Molecular Profiles:
#'      * molecularProfiles - Produce a molecular profiles dataset for a given
#'      study identifier ('studyId')
#'
#' @param projection character(default: "SUMMARY") Specify the projection
#'   type for data retrieval for details see API documentation
#'
#' @export
molecularProfiles <- function(api, studyId = NA_character_,
    projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    projection <- match.arg(projection)
    mols <- .invoke_fun(
        api = api, name = "getAllMolecularProfilesInStudyUsingGET",
        use_cache = TRUE, studyId = studyId, projection = projection
    )
    cmols <- httr::content(mols)
    if (projection %in% c("SUMMARY", "ID"))
        dplyr::bind_rows(cmols)
    else
        cmols
}

#' @name cBioPortal
#'
#' @section Molecular Profiles:
#'     * molecularData - Produce a dataset of molecular profile data based on
#'     `molecularProfileId`, `entrezGeneIds`, and `sampleIds`
#'
#' @export
molecularData <- function(api, molecularProfileId = NA_character_,
    entrezGeneIds = NULL, sampleIds = NULL, check = TRUE)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validMolProf <- .checkIdValidity(api, element = molecularProfileId,
        ename = "molecularProfileId", check = check)
    if (!validMolProf)
        stop("Provide a valid 'molecularProfileId' from 'molecularProfiles()'")
    if (is.null(entrezGeneIds))
        stop("Provide a character vector of 'entrezGeneIds'")
    if (is.null(sampleIds))
        stop("Provide a character vector of 'sampleIds'")

    mutation <- grepl("mutation", molecularProfileId)
    endpoint <- if (mutation) "fetchMutationsInMolecularProfileUsingPOST"
        else "fetchAllMolecularDataInMolecularProfileUsingPOST"

    byGene <- .invoke_bind(api,
        endpoint,
        use_cache = TRUE,
        molecularProfileId = molecularProfileId,
        entrezGeneIds = sort(entrezGeneIds),
        sampleIds = sort(sampleIds)
    )

    if ("message" %in% names(byGene) || !length(byGene)) {
        msg <- byGene[["message"]]
        if (length(msg))
            warning(msg)
        else
            warning(
                "No data found for molecularProfileId: ", molecularProfileId
            )
        dplyr::tibble()
    } else {
        byGene
    }
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * searchOps - Search through API operations with a keyword
#'
#' @export
searchOps <- function(api, keyword) {
    grep(keyword, names(AnVIL::operations(api)),
        value = TRUE, ignore.case = TRUE)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * geneTable - Get a table of all genes by 'entrezGeneId' or
#'     'hugoGeneSymbol'
#'
#' @param ... Additional arguments to lower level API functions
#'
#' @export
geneTable <- function(api, ...) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    gres <- .invoke_fun(api, "getAllGenesUsingGET", TRUE, ...)
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
samplesInSampleLists <-
    function(api, sampleListIds = NA_character_, check = TRUE) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validSLI <- .checkIdValidity(api, element = sampleListIds,
        ename = "sampleListId", use_cache = TRUE, check = check)
    if (!validSLI)
        stop("Provide valid 'sampleListIds' from 'sampleLists()'")

    sampleListIds <- sort(sampleListIds)
    sampleListIds <- stats::setNames(sampleListIds, sampleListIds)

    digi <- .inputDigest(match.call(), "samplesInSampleLists")
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
        meta <- structure(vector("list", length(sampleListIds)),
            .Names = sampleListIds)
        res <- lapply(sampleListIds, function(x) {
            res <- .invoke_fun(
                api, "getSampleListUsingGET", TRUE, sampleListId = x
            )
            res2 <- httr::content(res)
            meta[[x]] <<- res2[names(res2) != "sampleIds"]
            unlist(res2[["sampleIds"]])
        })
        res <- IRanges::CharacterList(res)
        meta <- dplyr::bind_rows(meta)
        metadata(res) <- meta
        save(res, file = cacheloc)
    }
    res
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * sampleLists - obtain all `sampleListIds` for a particular `studyId`
#'
#' @export
sampleLists <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    .invoke_bind(api, "getAllSampleListsInStudyUsingGET", TRUE,
        studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * allSamples - obtain all samples within a particular `studyId`
#'
#' @export
allSamples <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    .invoke_bind(api, "getAllSamplesInStudyUsingGET", TRUE, studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * genePanels - Show all available gene panels
#'
#' @export
genePanels <- function(api) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    .invoke_bind(api, "getAllGenePanelsUsingGET", TRUE)
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getGenePanels - Obtain the gene panel for a particular 'genePanelId'
#'
#' @export
getGenePanel <- function(api, genePanelId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    validGP <- .checkIdValidity(api, element = genePanelId,
        ename = "genePanelId")
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    res <- .invoke_fun(api, "getGenePanelUsingGET", TRUE,
        genePanelId = genePanelId)
    res <- httr::content(res)[["genes"]]
    dplyr::bind_rows(res)
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * genePanelMolecular - get gene panel data for a paricular
#'     `molecularProfileId` and `sampleListId` combination
#'
#' @param sampleListId character(1) A sample list identifier as obtained from
#'     `sampleLists()``
#'
#' @export
genePanelMolecular <-
    function(api, molecularProfileId = NA_character_, sampleListId = NULL,
        sampleIds = NULL)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    validMolProf <- .checkIdValidity(api,
        element = molecularProfileId, ename = "molecularProfileId")
    if (!validMolProf)
        stop("Provide a valid 'molecularProfileId' from 'molecularProfiles()'")

    if (!is.null(sampleListId))
        .invoke_bind(api, "getGenePanelDataUsingPOST", use_cache = TRUE,
            molecularProfileId = molecularProfileId,
            sampleListId = list(sampleListId = sampleListId)
        )
    else if (!is.null(sampleIds))
        .invoke_bind(api, "getGenePanelDataUsingPOST", use_cache = TRUE,
            molecularProfileId = molecularProfileId,
            sampleIds = list(sampleIds = sort(sampleIds))
        )
    else
        stop("Provide either 'sampleIds' or a 'sampleListId'")
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getGenePanelMolecular - get gene panel data for a combination of
#'     `molecularProfileId` and `sampleListId` vectors
#'
#' @export
getGenePanelMolecular <-
    function(api, molecularProfileIds = NA_character_, sampleIds)
{
    validMolProf <- .checkIdValidity(api,
        element = molecularProfileIds, ename = "molecularProfileId")
    if (!validMolProf)
        stop(
            paste0("Provide multiple valid 'molecularProfileIds' from",
            " 'molecularProfiles()'")
        )
    if (missing(sampleIds))
        stop(
            paste0("Provide valid 'sampleIds' from 'samplesInSampleLists()'",
            " or 'allSamples()'")
        )

    SampMolIds <- S4Vectors::expand.grid(
        molecularProfileId = sort(molecularProfileIds),
        sampleId = sort(sampleIds)
    )
    SampMolIds <- SampMolIds[order(SampMolIds[["molecularProfileId"]]), ]

    .invoke_bind(
        api = api,
        name = "fetchGenePanelDataInMultipleMolecularProfilesUsingPOST",
        use_cache = TRUE,
        sampleMolecularIdentifiers = SampMolIds
    )
}

#' @name cBioPortal
#'
#' @section Sample Data:
#'     * getSampleInfo - Obtain sample metadata for a particular `studyId` or
#'     `sampleListId`
#' @export
getSampleInfo <-
    function(api, studyId = NA_character_, sampleListIds = NULL,
        projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    projection <- match.arg(projection)
    if (!is.null(sampleListIds))
        queryobj <- list(sampleListIds = sampleListIds)
    else
        queryobj <- list(sampleIdentifiers =
            as.data.frame(
                allSamples(api, studyId)[, c("sampleId", "studyId")]
            )
        )

    .invoke_bind(api = api, name = "fetchSamplesUsingPOST", use_cache = TRUE,
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
#' @param check logical(1) Whether to check the inputs against values from the
#'     API (i.e., for 'studyId', 'genePanelId', 'molecularProfileId', and
#'     'sampleListId')
#'
#' @examples
#'
#' getDataByGenePanel(cbio, studyId = "acc_tcga", genePanelId = "IMPACT341",
#'    molecularProfileId = "acc_tcga_rppa")
#'
#' @export
getDataByGenePanel <-
    function(api, studyId = NA_character_, genePanelId = NA_character_,
        molecularProfileId = NULL, sampleListId = NULL, check = TRUE)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId",
        check = check)
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    validGP <- .checkIdValidity(api, element = genePanelId,
        ename = "genePanelId", check = check)
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    if (!is.null(sampleListId))
        samples <- samplesInSampleLists(api, sampleListId, check = check)[[1L]]
    else
        samples <- allSamples(api, studyId)[["sampleId"]]

    panel <- getGenePanel(api, genePanelId = genePanelId)
    molecularData <- molecularData(api = api,
        molecularProfileId = molecularProfileId,
        entrezGeneIds = panel[["entrezGeneId"]],
        sampleIds = samples, check = check)

    suppressMessages({
        dplyr::left_join(molecularData, panel)
    })
}
