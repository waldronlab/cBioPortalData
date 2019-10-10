utils::globalVariables(c("clinicalAttributeId", "value", "sampleId"))

.invoke_fun <- function(api, name, do_cache = FALSE, ...) {
    if (!is(api, "cBioPortal"))
        stop("Provide a 'cBioPortal' class API object")
    ops <- names(AnVIL::operations(api))
    if (!name %in% ops)
        stop("<internal> operation name not found in API")

    if (do_cache) {
        .dollarCache(list(api, name), ...)
    } else {
        do.call(`$`, list(api, name))(...)
    }
}

.dollarCache <- function(appname, ...) {
    if (!is.list(appname))
        stop("<internal> Provide a list input as 'api$name'")
    digi <- digest::digest(list(appname, ...))
    loc <- .getHashCache(digi)
    if (file.exists(loc)) {
        load(loc)
    } else {
        op <- do.call(`$`, appname)(...)
        save(op, file = loc)
    }
    op
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_bind <- function(api, name, do_cache = FALSE, ...) {
    .bind_content(.invoke_fun(api, name, do_cache, ...))
}

#' @name cBioPortal-class
#'
#' @title A class for representing the cBioPortal API
#'
#' @description The cBioPortal class is a product of the cBioPortal API that
#'     inherits from the `Service` class
#'
#' @importFrom methods new
#'
#' @seealso AnVIL::Service
#'
#' @examples
#'
#' cbio <- cBioPortal()
#'
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
#'
#' getStudies(cbio = cc)
#'
#' searchOps(cbio = cc, keyword = "molecular")
#'
#' molecularProfiles(cbio = cc, studyId = "acc_tcga")
#'
#' molecularSlice(cbio = cc, molecularProfileId = "acc_tcga_rna_seq_v2_mrna",
#'     entrezGeneIds = c(1, 2),
#'     sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#' )
#'
#' sampleLists(cbio = cc, studyId = "acc_tcga")
#'
#' samplesInSampleLists(cbio = cc,
#'     sampleListIds = c("acc_tcga_rppa", "acc_tcga_cnaseq")
#' )
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
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    digi <- .inputDigest(match.call())
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
        query <- .invoke_fun(cbio, "getAllStudiesUsingGET")
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
#' @section clinicalData:
#'      Obtain clinical data for a particular study identifier
#'
#' @export
clinicalData <- function(cbio, studyId = NA_character_) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    digcall <- match.call()
    where <- function(name, env = parent.frame()) {
        if (identical(env, emptyenv()))
            stop("Can't find ", name, call. = FALSE)
        else if (exists(name, envir = env, inherits = FALSE))
            env
        else
            where(name, parent.env(env))
    }

    if (sys.nframe() > 1L) {
        a <- as.list(digcall)
        b <- as.list(match.call(sys.function(sys.parent(1L)),
            call = sys.call(1L), envir = where('cbio')))
        comm <- Filter(nchar, intersect(names(a), names(b)))
        a[comm] <- b[comm]
        digcall <- as.call(a)
    }

    digi <- .inputDigest(digcall)
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
        pttable <- .invoke_bind(cbio,
            "getAllPatientsInStudyUsingGET", studyId = studyId)
        ptrow <- lapply(pttable[["patientId"]], function(pt) {
            .invoke_bind(cbio,
                "getAllClinicalDataOfPatientInStudyUsingGET",
                studyId = studyId, patientId = pt)
        })
        save(ptrow, file = cacheloc)
    }
    clin <- dplyr::bind_rows(ptrow)
    tidyr::spread(clin, clinicalAttributeId, value)
}

#' @name cBioPortal
#'
#' @section Molecular Profiles:
#'      * molecularProfiles - Produce a molecular profiles dataset for a given
#'      study identifier
#'
#' @param projection character(default: "SUMMARY") Specify the projection
#'   type for data retrieval for details see API documentation
#'
#' @inheritParams cBioPortal
#'
#' @export
molecularProfiles <- function(cbio, studyId = NA_character_,
    projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    projection <- match.arg(projection)
    mols <- .invoke_fun(cbio, "getAllMolecularProfilesInStudyUsingGET",
        do_cache = TRUE, studyId = studyId, projection = projection)
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
molecularSlice <- function(cbio, molecularProfileId = NA_character_,
    entrezGeneIds = NULL, sampleIds = NULL, check = TRUE)
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validMolProf <- .checkIdValidity(cbio, element = molecularProfileId,
        ename = "molecularProfileId", check = check)
    if (!validMolProf)
        stop("Provide a valid 'molecularProfileId' from 'molecularProfiles()'")
    if (is.null(entrezGeneIds))
        stop("Provide a character vector of 'entrezGeneIds'")
    if (is.null(sampleIds))
        stop("Provide a character vector of 'sampleIds'")

    byGene <- .invoke_bind(cbio,
        "fetchAllMolecularDataInMolecularProfileUsingPOST",
        do_cache = TRUE,
        molecularProfileId = molecularProfileId,
        entrezGeneIds = sort(entrezGeneIds),
        sampleIds = sort(sampleIds)
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
    grep(keyword, names(AnVIL::operations(cbio)),
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
geneTable <- function(cbio, ...) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    gres <- .invoke_fun(cbio, "getAllGenesUsingGET", TRUE, ...)
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
    function(cbio, sampleListIds = NA_character_, check = TRUE) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validSLI <- .checkIdValidity(cbio, element = sampleListIds,
        ename = "sampleListIds", do_cache = TRUE, check = check)
    if (!validSLI)
        stop("Provide valid 'sampleListIds' from 'sampleLists()'")

    sampleListIds <- sort(sampleListIds)
    sampleListIds <- stats::setNames(sampleListIds, sampleListIds)

    digi <- .inputDigest(match.call())
    cacheloc <- .getHashCache(digi)
    if (file.exists(cacheloc)) {
        load(cacheloc)
    } else {
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
        save(res, file = cacheloc)
    }
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
sampleLists <- function(cbio, studyId = NA_character_) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    .invoke_bind(cbio, "getAllSampleListsInStudyUsingGET", TRUE,
        studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * allSamples - obtain all samples within a particular `studyId`
#'
#' @export
allSamples <- function(cbio, studyId = NA_character_) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    .invoke_bind(cbio, "getAllSamplesInStudyUsingGET", TRUE, studyId = studyId)
}

#' @name cBioPortal
#'
#' @section API Metadata:
#'     * genePanels - Show all available gene panels
#'
#' @export
genePanels <- function(cbio) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    .invoke_bind(cbio, "getAllGenePanelsUsingGET", TRUE)
}

#' @name cBioPortal
#'
#' @section Gene Panels:
#'     * getGenePanels - Obtain the gene panel for a particular 'genePanelId'
#'
#' @export
getGenePanel <- function(cbio, genePanelId = NA_character_) {
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    validGP <- .checkIdValidity(cbio, element = genePanelId,
        ename = "genePanelId")
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    res <- .invoke_fun(cbio, "getGenePanelUsingGET", TRUE,
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
    function(cbio, molecularProfileId = NA_character_, sampleListId = NULL,
        sampleIds = NULL)
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    validMolProf <- .checkIdValidity(cbio,
        element = molecularProfileId, ename = "molecularProfileId")
    if (!validMolProf)
        stop("Provide a valid 'molecularProfileId' from 'molecularProfiles()'")

    if (!is.null(sampleListId))
        .invoke_bind(cbio, "getGenePanelDataUsingPOST", cache = TRUE,
            molecularProfileId = molecularProfileId,
            sampleListId = list(sampleListId = sampleListId)
        )
    else if (!is.null(sampleIds))
        .invoke_bind(cbio, "getGenePanelDataUsingPOST", cache = TRUE,
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
    function(cbio, molecularProfileIds = NA_character_, sampleIds)
{
    validMolProf <- .checkIdValidity(cbio,
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
    .invoke_bind(cbio,
        "fetchGenePanelDataInMultipleMolecularProfilesUsingPOST",
        TRUE,
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
    function(cbio, studyId = NA_character_, sampleListIds = NULL,
        projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    projection <- match.arg(projection)
    if (!is.null(sampleListIds))
        queryobj <- list(sampleListIds = sampleListIds)
    else
        queryobj <- list(sampleIdentifiers =
            as.data.frame(
                allSamples(cbio, studyId)[, c("sampleId", "studyId")]
            )
        )

    .invoke_bind(cbio, "fetchSamplesUsingPOST", TRUE,
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
#' @param by character(1) Either 'entrezGeneId' or 'hugoGeneSymbol' for row
#'     metadata
#'
#' @param check logical(1) Whether to check the inputs against values from the
#'     API (i.e., for 'studyId', 'genePanelId', 'molecularProfileId', and
#'     'sampleListId')
#'
#' @export
getDataByGenePanel <-
    function(cbio, studyId = NA_character_, genePanelId = NA_character_,
        by = c("entrezGeneId", "hugoGeneSymbol"),
        molecularProfileId = NULL, sampleListId = NULL, check = TRUE)
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")
    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId",
        check = check)
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    validGP <- .checkIdValidity(cbio, element = genePanelId,
        ename = "genePanelId", check = check)
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    by <- match.arg(by)
    if (!is.null(sampleListId))
        samples <- samplesInSampleLists(cbio, sampleListId, check = check)[[1L]]
    else
        samples <- allSamples(cbio, studyId)[["sampleId"]]

    panel <- getGenePanel(cbio, genePanelId = genePanelId)
    molecularData <- molecularSlice(cbio = cbio,
        molecularProfileId = molecularProfileId,
        entrezGeneIds = panel[["entrezGeneId"]],
        sampleIds = samples, check = check)

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
    function(cbio, by, genePanelId, studyId, molecularProfileIds, sampleListId,
        check)
{
    if (is.null(molecularProfileIds)) {
        molecularProfileIds <-
            molecularProfiles(cbio, studyId)[["molecularProfileId"]]
    } else { check <- TRUE }

    molecularProfileIds <- stats::setNames(molecularProfileIds,
        molecularProfileIds)

    expers <- lapply(molecularProfileIds, function(molprof) {
        moldata <- getDataByGenePanel(cbio, by = by,
            genePanelId = genePanelId, studyId = studyId,
            molecularProfileId = molprof, sampleListId = sampleListId,
            check = check)
        moldata <- as.data.frame(moldata)
        rownames(moldata) <- moldata[[by]]
        moldata <- data.matrix(moldata[, names(moldata) != by])
        SummarizedExperiment(moldata)
    })
    as(Filter(length, expers), "List")
}

std.args <- function(call, formals) {
    callargs <- as.list(call)[-1]
    toadd <- setdiff(names(formals), names(call))
    call[toadd] <- formals[toadd]
    call
}

match.args <- function(fun, call, ...) {
    funfor <- formals(fun)
    exargs <- intersect(names(funfor), names(call))
    c(as.list(call)[-1][exargs], ...)
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
#' @param idConvert function(default: `identity()`) A function to process
#'     identifiers for matching patients to samples. It catches TCGA samples
#'     automatically and uses `TCGAutils::TCGAbarcode` instead.
#'
#' @examples
#'
#' cb <- cBioPortal()
#' cBioPortalData(cb, by = "hugoGeneSymbol", studyId = "acc_tcga",
#'     genePanelId = "IMPACT341",
#'     molecularProfileIds = c("acc_tcga_rppa", "acc_tcga_linear_CNA")
#' )
#'
#' @export
cBioPortalData <-
    function(cbio, studyId = NA_character_,
        genePanelId = NA_character_,
        molecularProfileIds = NULL,
        sampleListId = NULL,
        by = c("entrezGeneId", "hugoGeneSymbol"),
        idConvert = identity
    )
{
    if (missing(cbio))
        stop("Provide a valid 'cbio' from 'cBioPortal()'")

    validStudy <- .checkIdValidity(cbio, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    validGP <-
        .checkIdValidity(cbio, element = genePanelId, ename = "genePanelId")
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    by <- match.arg(by)

    call <- std.args(match.call(), formals())
    exargs <- match.args(.portalExperiments, call, check = FALSE)
    explist <- do.call(.portalExperiments, exargs)

    explist <- as(explist, "ExperimentList")

    clin <- do.call(clinicalData, match.args(clinicalData, call))
    clin <- as.data.frame(clin)
    rownames(clin) <- clin[["patientId"]]
    if (all(startsWith(rownames(clin), "TCGA")))
        idConvert <- TCGAutils::TCGAbarcode

    sampmap <- try({
        TCGAutils::generateMap(experiments = explist,
            colData = clin, idConverter = idConvert)
    }, silent = TRUE)

    if (is(sampmap, "try-error"))
        idConvert <- .generateIdConvert(
            unlist(colnames(explist), use.names = FALSE),
            rownames(clin)
        )

    sampmap <- TCGAutils::generateMap(experiments = explist,
        colData = clin, idConverter = idConvert)

    MultiAssayExperiment(explist, clin, sampmap)
}
