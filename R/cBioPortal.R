.parse_token <- function(token_file) {
    token <- try({
        as.character(read.dcf(token_file, fields = "token"))
    })
    if (is(token, "try-error"))
        token <- readLines(token_file)
    token
}

.handle_token <- function(token) {
    if (file.exists(token))
        token <- .parse_token(token)
    else if (grepl(.Platform$file.sep, token, fixed = TRUE))
        stop("The token filepath is not valid")
    token <- gsub("token: ", "", token)
    c(Authorization = paste("Bearer", token))
}

#' The R interface to the cBioPortal API Data Service
#'
#' @description This section of the documentation lists the functions that allow
#'   users to access the cBioPortal API. The main representation of the API can
#'   be obtained from the `cBioPortal` function. The supporting functions listed
#'   here give access to specific parts of the API and allow the user to explore
#'   the API with individual calls. Many of the functions here are listed for
#'   documentation purposes and are recommended for advanced usage only. Users
#'   should only need to use the `cBioPortalData` main function to obtain data.
#'
#' @param api An API object of class `cBioPortal` from the `cBioPortal` function
#'
#' @param hostname `character(1)` The internet location of the service (default:
#'   'www.cbioportal.org')
#'
#' @param protocol `character(1)` The internet protocol used to access the
#'   hostname (default: 'https')
#'
#' @param api. `character(1)` The directory location of the API protocol within
#'   the hostname (default: '/api/v2/api-docs')
#'
#' @param token `character(1)` The Authorization Bearer token e.g.,
#'   "63eba81c-2591-4e15-9d1c-fb6e8e51e35d" or a path to text file.
#'
#' @param studyId `character(1)` Indicates the "studyId" as taken from
#'   `getStudies`
#'
#' @param buildReport `logical(1)` Indicates whether to append the build
#'   information to the `getStudies()` table (default FALSE)
#'
#' @param keyword `character(1)` Keyword or pattern for searching through
#'   available operations
#'
#' @param molecularProfileId `character(1)` Indicates a molecular profile ID
#'
#' @param molecularProfileIds `character()` A vector of molecular profile IDs
#'
#' @param entrezGeneIds `numeric()` A vector indicating entrez gene IDs
#'
#' @param sampleIds `character()` Sample identifiers
#'
#' @param genes `character()` Either Entrez gene identifiers or Hugo gene
#'   symbols. When included, the 'by' argument indicates the type of identifier
#'   provided and 'genePanelId' is ignored. Preference is given to Entrez IDs
#'   due to faster query responses.
#'
#' @param genePanelId `character(1)` Identifies the gene panel, as obtained from
#'   the `genePanels` function
#'
#' @param by `character(1)` Either 'entrezGeneId' or 'hugoGeneSymbol' for row
#'   metadata (default: 'entrezGeneId')
#'
#' @return
#'   * cBioPortal: An API object of class 'cBioPortal'
#'   * cBioPortalData: A data object of class 'MultiAssayExperiment'
#'
#' @importFrom AnVIL Service
#'
#' @examples
#' cbio <- cBioPortal()
#'
#' @export
cBioPortal <- function(
    hostname = "www.cbioportal.org",
    protocol = "https",
    api. = "/api/v2/api-docs",
    token = character()
) {
    if (length(token))
        token <- .handle_token(token)

    apiUrl <- paste0(protocol, "://", hostname, api.)
    service <- withCallingHandlers({
        Service(
            service = "cBioPortal",
            host = hostname,
            config = httr::config(
                ssl_verifypeer = 0L, ssl_verifyhost = 0L, http_version = 0L
            ),
            authenticate = FALSE,
            api_reference_url = apiUrl,
            api_reference_md5sum = "7314de5c5e8056e4e07b411b3e5a0cb9",
            api_reference_headers = token,
            package = "cBioPortalData",
            schemes = protocol
        )
    }, warning = function(w) {
        if (!grepl("incomplete final line", w))
            warning(w)
        invokeRestart("muffleWarning")
    })
    .cBioPortal(
        service, api_header = token
    )
}

#' @importFrom BiocBaseUtils checkInstalled
.loadReportData <- function() {
    checkInstalled("jsonlite")
    denv <- new.env(parent = emptyenv())
    api_build <- system.file(
        "extdata", "api", "api_build.json",
        package = "cBioPortalData", mustWork = TRUE
    ) |>
        jsonlite::fromJSON() |>
        as.data.frame()
    pack_build <- system.file(
        "extdata", "pack", "pack_build.json",
        package = "cBioPortalData", mustWork = TRUE
    ) |>
        jsonlite::fromJSON() |>
        as.data.frame()
    denv[["api_build"]] <- api_build
    denv[["pack_build"]] <- pack_build

    denv
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * getStudies: Obtain a table of studies and associated metadata and
#'   optionally include a `buildReport` status (default FALSE) for each
#'   study. When enabled, the 'api_build' and 'pack_build' columns will
#'   be added to the table and will show if `MultiAssayExperiment` objects
#'   can be generated for that particular study identifier (`studyId`). The
#'   'api_build' column corresponds to datasets obtained with
#'   `cBioPortalData` and the 'pack_build' column corresponds to datsets
#'   loaded via `cBioDataPack`.
#'
#' @examples
#' getStudies(api = cbio)
#'
#' @export
getStudies <- function(api, buildReport = FALSE) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    query <- .invoke_fun(api, "getAllStudiesUsingGET")
    studies <- httr::content(query)
    studies <- lapply(studies, function(x) {
        if (is.null(x[["pmid"]]))
            x[["pmid"]] <- NA_character_
        if (is.null(x[["citation"]]))
            x[["citation"]] <- NA_character_
        x
    })
    studytable <- dplyr::bind_rows(studies)

    if (buildReport) {
        denv <- .loadReportData()
        suppressMessages({
            studytable <- dplyr::left_join(studytable, denv[["api_build"]])
            studytable <- dplyr::left_join(studytable, denv[["pack_build"]])
        })
    }

    studytable
}

#' @rdname cBioPortal
#'
#' @section Patient Data:
#' * clinicalData - Obtain clinical data for a particular study identifier
#'   ('studyId')
#'
#' @examples
#' clinicalData(cbio, "acc_tcga")
#'
#' @export
clinicalData <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    atts <- .invoke_bind(
        api = api, name = "fetchAllClinicalDataInStudyUsingPOST",
        clinicalDataType = "PATIENT", studyId = studyId
    )
    att_tab <- tidyr::pivot_wider(data = atts, id_cols = "patientId",
        names_from = "clinicalAttributeId", values_from = "value")
    clin <- .invoke_bind(
        api = api, name = "getAllClinicalDataInStudyUsingGET",
        studyId = studyId
    )
    clin_tab <- tidyr::pivot_wider(data = clin,
        id_cols = c("patientId", "sampleId"),
        names_from = "clinicalAttributeId", values_from = "value")
    dplyr::full_join(att_tab, clin_tab, by = "patientId")
}

#' @rdname cBioPortal
#'
#' @section Molecular Profiles:
#' * molecularProfiles - Produce a molecular profiles dataset for a given
#'   study identifier ('studyId')
#'
#' @param projection `character(1)` (default: "SUMMARY") Specify the projection
#'   type for data retrieval for details see API documentation
#'
#' @examples
#' molecularProfiles(cbio, "acc_tcga")
#'
#' @export
molecularProfiles <- function(api, studyId = NA_character_,
    projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    projection <- match.arg(projection)
    mols <- .invoke_fun(
        api = api, name = "getAllMolecularProfilesInStudyUsingGET",
        studyId = studyId, projection = projection
    )
    cmols <- httr::content(mols)
    if (projection %in% c("SUMMARY", "ID"))
        dplyr::bind_rows(cmols)
    else
        cmols
}

.sampleMolIds <- function(molecularProfileIds, sampleIds)
{
    SampMolIds <- S4Vectors::expand.grid(
        molecularProfileId = sort(molecularProfileIds),
        sampleId = sort(sampleIds)
    )
    SampMolIds[order(SampMolIds[["molecularProfileId"]]), ]
}

.FilterLengthWarn <- function(datalist) {
    datanames <- names(datalist)
    for (dname in datanames) {
        byG <- datalist[[dname]]
        if ("message" %in% names(byG) || !length(byG)) {
            msg <- byG[["message"]]
            if (length(msg)) {
                warning(dname, ": ", msg, call. = FALSE)
                datalist[[dname]] <- dplyr::tibble()
            }
        }
    }
    empty <- vapply(datalist, function(adat) !length(adat), logical(1L))
    if (any(empty)) {
        enames <- paste(datanames[empty], collapse = ", ")
        warning("No data found for molecularProfileId: ", enames, call. = FALSE)
    }
    ## remove empty responses (e.g., in ov_tcga_pub_mirna)
    Filter(length, datalist)
}

#' @rdname cBioPortal
#'
#' @section Molecular Data:
#' * fetchData - A convenience function to download both mutation and
#'   molecular data with `molecularProfileId`, `entrezGeneIds`, and
#'   `sampleIds`
#'
#' @examples
#' fetchData(
#'     api = cbio, studyId = "acc_tcga",
#'     molecularProfileIds = c(
#'         "acc_tcga_mutations", "acc_tcga_gistic", "acc_tcga_rppa"
#'     ),
#'     entrezGeneIds = 1:1000,
#'     sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#' )
#' @export
fetchData <-
    function(
        api, studyId, molecularProfileIds = NA_character_,
        entrezGeneIds = NULL, sampleIds = NULL
) {
    if (missing(studyId))
        stop("Provide a valid 'studyId' from 'cBioPortal()'")
    byGeneList <- vector("list", length(molecularProfileIds))
    names(byGeneList) <- molecularProfileIds

    all_profs <- molecularProfiles(api = api, studyId = studyId)
    mut_profs <- all_profs[
        all_profs[["molecularAlterationType"]] == "MUTATION_EXTENDED",
    ] |>
        dplyr::pull("molecularProfileId")
    mutation <- molecularProfileIds %in% mut_profs
    mutationList <- mutationData(
        api = api, molecularProfileIds = molecularProfileIds[mutation],
        entrezGeneIds = entrezGeneIds, sampleIds = sampleIds
    )

    molecularProfileIds <- molecularProfileIds[!mutation]
    dcn_molprofs <- all_profs[
        all_profs[["molecularAlterationType"]] == "COPY_NUMBER_ALTERATION" &
            all_profs[["datatype"]] == "DISCRETE",
    ] |>
        dplyr::pull("molecularProfileId")
    dcn <- molecularProfileIds %in% dcn_molprofs

    dcnList <- copyNumberData(
        api = api, molecularProfileIds = molecularProfileIds[dcn],
        entrezGeneIds = entrezGeneIds, sampleIds = sampleIds
    )

    molecularList <- molecularData(
        api = api, molecularProfileIds = molecularProfileIds[!dcn],
        entrezGeneIds = entrezGeneIds, sampleIds = sampleIds
    )
    byGeneList <- c(mutationList, dcnList, molecularList)
    .FilterLengthWarn(byGeneList)
}

#' @rdname cBioPortal
#'
#' @section Molecular Data:
#' * mutationData - Produce a dataset of mutation data using
#'   `molecularProfileId`, `entrezGeneIds`, and `sampleIds`
#'
#' @examples
#' mutationData(
#'     api = cbio,
#'     molecularProfileIds = "acc_tcga_mutations",
#'     entrezGeneIds = 1:1000,
#'     sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#' )
#'
#' @export
mutationData <- function(api, molecularProfileIds = NA_character_,
    entrezGeneIds = NULL, sampleIds = NULL)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    if (is.null(entrezGeneIds))
        stop("Provide a character vector of 'entrezGeneIds'")
    if (is.null(sampleIds))
        stop("Provide a character vector of 'sampleIds'")

    if (length(molecularProfileIds) > 1L) {
        SampMolIds <- .sampleMolIds(molecularProfileIds, sampleIds)
    }

    if (length(molecularProfileIds) == 1L) {
        endpoint <- "fetchMutationsInMolecularProfileUsingPOST"
        byGene <- .invoke_bind(
            api, endpoint,
            molecularProfileId = molecularProfileIds,
            entrezGeneIds = sort(entrezGeneIds),
            sampleIds = sort(sampleIds)
        )
    } else if (length(molecularProfileIds)) {
        endpoint <- "fetchMutationsInMultipleMolecularProfilesUsingPOST"
        byGene <- .invoke_bind(
            api, endpoint,
            molecularProfileIds = molecularProfileIds,
            sampleMolecularIdentifiers = SampMolIds
        )
    }
    if (!length(molecularProfileIds) || !length(byGene))
        structure(
            vector("list", length(molecularProfileIds)),
            .Names = molecularProfileIds
        )
    else
        split(byGene, byGene[["molecularProfileId"]])
}

#' @rdname cBioPortal
#'
#' @section Molecular Data:
#' * molecularData - Produce a dataset of molecular profile data based on
#'   `molecularProfileId`, `entrezGeneIds`, and `sampleIds`
#'
#' @examples
#' molecularData(
#'     api = cbio,
#'     molecularProfileIds = c("acc_tcga_rna_seq_v2_mrna", "acc_tcga_rppa"),
#'     entrezGeneIds = 1:100,
#'     sampleIds = c("TCGA-OR-A5J1-01", "TCGA-OR-A5J2-01")
#' )
#'
#' @export
molecularData <- function(api, molecularProfileIds = NA_character_,
    entrezGeneIds = NULL, sampleIds = NULL)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    if (is.null(entrezGeneIds))
        stop("Provide a character vector of 'entrezGeneIds'")
    if (is.null(sampleIds))
        stop("Provide a character vector of 'sampleIds'")

    if (length(molecularProfileIds) == 1L) {
        endpoint <- "fetchAllMolecularDataInMolecularProfileUsingPOST"
        byGene <- .invoke_bind(api,
            endpoint,
            molecularProfileId = molecularProfileIds,
            entrezGeneIds = sort(entrezGeneIds),
            sampleIds = sort(sampleIds)
        )
    } else if (length(molecularProfileIds)) {
        endpoint <- "fetchMolecularDataInMultipleMolecularProfilesUsingPOST"
        byGene <- .invoke_bind(api,
            endpoint,
            projection = "SUMMARY",
            entrezGeneIds = sort(entrezGeneIds),
            sampleMolecularIdentifiers = .sampleMolIds(
                molecularProfileIds, sampleIds
            )
        )
    }
    if (!length(molecularProfileIds) || !length(byGene))
        structure(
            vector("list", length(molecularProfileIds)),
            .Names = molecularProfileIds
        )
    else
        split(byGene, byGene[["molecularProfileId"]])
}

#' @rdname cBioPortal
#'
#' @section Copy Number Data:
#' * copyNumberData - Produce a dataset of copy number data based on
#'   `molecularProfileId`, `sampleListId`, `discreteCopyNumberEventType`, and
#'   `projection`
#'
#' @param discreteCopyNumberEventType `character(1)` The copy number event type
#'   to filter on. Must be one of "HOMDEL_AND_AMP" (default), "HOMDEL", "AMP",
#'   "GAIN", "HETLOSS", "DIPLOID", or "ALL"
#'
#' @examples
#' ## obtain molecularProfileId for discrete copy number alteration data
#' molecularProfiles(cbio, "acc_tcga") |>
#'     dplyr::filter(
#'         molecularAlterationType == "COPY_NUMBER_ALTERATION" &
#'         datatype == "DISCRETE"
#'     )
#'
#' copyNumberData(
#'     api = cbio,
#'     molecularProfileIds = "acc_tcga_gistic",
#'     entrezGeneIds = 25,
#'     sampleListId = "acc_tcga_all"
#' )
#'
#' @export
copyNumberData <- function(
    api, molecularProfileIds = NA_character_,
    entrezGeneIds = NULL,
    sampleIds = NULL, sampleListId = NULL,
    discreteCopyNumberEventType = c(
        "HOMDEL_AND_AMP", "HOMDEL", "AMP", "GAIN", "HETLOSS", "DIPLOID", "ALL"
    ),
    projection = c("SUMMARY", "ID", "DETAILED", "META")
) {
    discreteCopyNumberEventType <- match.arg(discreteCopyNumberEventType)
    projection <- match.arg(projection)
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    if (is.null(entrezGeneIds))
        stop("Provide a character vector of 'entrezGeneIds'")
    if (is.null(sampleListId) && is.null(sampleIds))
        stop("Provide either a 'sampleListId' or 'sampleIds'")

    if (!length(molecularProfileIds) || all(is.na(molecularProfileIds)))
        return(
            structure(
                vector("list", length(molecularProfileIds)),
                .Names = molecularProfileIds
            )
        )

    names(molecularProfileIds) <- molecularProfileIds
    lapply(
        molecularProfileIds,
        function(molecularProfileId) {
            endpoint <- "fetchDiscreteCopyNumbersInMolecularProfileUsingPOST"
            .invoke_bind(
                api, endpoint,
                molecularProfileId = molecularProfileId,
                entrezGeneIds = sort(entrezGeneIds),
                sampleListId = sampleListId,
                sampleIds = sort(sampleIds),
                discreteCopyNumberEventType = discreteCopyNumberEventType,
                projection = projection
            )
        }
    )
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * searchOps - Search through API operations with a keyword
#'
#' @examples
#' searchOps(api = cbio, keyword = "molecular")
#'
#' @export
searchOps <- function(api, keyword) {
    grep(keyword, names(AnVIL::operations(api)),
        value = TRUE, ignore.case = TRUE)
}

#' @rdname cBioPortal
#'
#' @section Sample Data:
#' * samplesInSampleLists - get all samples associated with a 'sampleListId'
#'
#' @param sampleListIds `character()` A vector of 'sampleListId' as obtained
#'   from `sampleLists`
#'
#' @examples
#' samplesInSampleLists(
#'     api = cbio,
#'     sampleListIds = c("acc_tcga_rppa", "acc_tcga_cnaseq")
#' )
#'
#' @export
samplesInSampleLists <-
    function(api, sampleListIds = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    sampleListIds <- sort(sampleListIds)
    sampleListIds <- setNames(sampleListIds, sampleListIds)

    meta <- structure(vector("list", length(sampleListIds)),
        .Names = sampleListIds)
    res <- lapply(sampleListIds, function(x) {
        res <- .invoke_fun(
            api, "getSampleListUsingGET", sampleListId = x
        )
        res2 <- httr::content(res)
        meta[[x]] <<- res2[names(res2) != "sampleIds"]
        unlist(res2[["sampleIds"]])
    })
    res <- IRanges::CharacterList(res)
    meta <- dplyr::bind_rows(meta)
    metadata(res) <- meta

    res
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * sampleLists - obtain all `sampleListIds` for a particular `studyId`
#'
#' @examples
#' sampleLists(api = cbio, studyId = "acc_tcga")
#'
#' @export
sampleLists <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    .invoke_bind(api, "getAllSampleListsInStudyUsingGET",
        studyId = studyId)
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * allSamples - obtain all samples within a particular `studyId`
#'
#' @export
allSamples <- function(api, studyId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    .invoke_bind(api, "getAllSamplesInStudyUsingGET",
        studyId = list(studyId = studyId))
}

#' @rdname cBioPortal
#'
#' @section Sample Data:
#' * getSampleInfo - Obtain sample metadata for a particular `studyId` or
#'   `sampleListId`
#'
#' @export
getSampleInfo <-
    function(api, studyId = NA_character_, sampleListIds = NULL,
        projection = c("SUMMARY", "ID", "DETAILED", "META"))
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    projection <- match.arg(projection)
    if (!is.null(sampleListIds))
        queryobj <- list(sampleListIds = sampleListIds)
    else
        queryobj <- list(sampleIdentifiers =
            as.data.frame(
                allSamples(api, studyId)[, c("sampleId", "studyId")]
            )
        )

    .invoke_bind(api = api, name = "fetchSamplesUsingPOST",
        projection = projection, sampleIdentifiers = queryobj
    )
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * genePanels - Show all available gene panels
#'
#' @examples
#' genePanels(cbio)
#'
#' @export
genePanels <- function(api) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    .invoke_bind(api, "getAllGenePanelsUsingGET")
}

#' @rdname cBioPortal
#'
#' @section Gene Panels:
#' * getGenePanels - Obtain the gene panel for a particular 'genePanelId'
#'
#' @examples
#' getGenePanel(cbio, "AmpliSeq")
#'
#' @export
getGenePanel <- function(api, genePanelId = NA_character_) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    res <- .invoke_fun(api, "getGenePanelUsingGET", genePanelId = genePanelId)
    res <- httr::content(res)[["genes"]]
    dplyr::bind_rows(res)
}

#' @rdname cBioPortal
#'
#' @section Gene Panels:
#' * genePanelMolecular - get gene panel data for a particular
#'   `molecularProfileId` and either a vector of `sampleListId` or `sampleId`
#'
#' @param sampleListId `character(1)` A sample list identifier as obtained from
#'     `sampleLists()`
#'
#' @export
genePanelMolecular <- function(
    api, molecularProfileId = NA_character_,
    sampleListId = NULL, sampleIds = NULL
) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    if (!is.null(sampleListId))
        .invoke_bind(api, "getGenePanelDataUsingPOST",
            molecularProfileId = molecularProfileId,
            sampleListId = list(sampleListId = sampleListId)
        )
    else if (!is.null(sampleIds))
        .invoke_bind(api, "getGenePanelDataUsingPOST",
            molecularProfileId = molecularProfileId,
            sampleIds = list(sampleIds = sort(sampleIds))
        )
    else
        stop("Provide either 'sampleIds' or a 'sampleListId'")
}

#' @rdname cBioPortal
#'
#' @section Gene Panels:
#' * getGenePanelMolecular - get gene panel data for multiple
#'   `molecularProfileId`s and a vector of `sampleIds`
#'
#' @export
getGenePanelMolecular <-
    function(api, molecularProfileIds = NA_character_, sampleIds)
{
    if (missing(sampleIds))
        stop(
            "Provide valid 'sampleIds' from 'samplesInSampleLists()'",
            " or 'allSamples()'"
        )

    SampMolIds <- S4Vectors::expand.grid(
        molecularProfileId = sort(molecularProfileIds),
        sampleId = sort(sampleIds)
    )
    SampMolIds <- SampMolIds[order(SampMolIds[["molecularProfileId"]]), ]

    .invoke_bind(
        api = api,
        name = "fetchGenePanelDataInMultipleMolecularProfilesUsingPOST",
        sampleMolecularIdentifiers =
            list(sampleMolecularIdentifiers = SampMolIds)
    )
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * geneTable - Get a table of all genes by 'entrezGeneId' and
#'   'hugoGeneSymbol'
#'
#' @param pageSize `numeric(1)` The number of rows in the table to return
#'
#' @param pageNumber `numeric(1)` The pagination page number
#'
#' @param ... Additional arguments to lower level API functions
#'
#' @export
geneTable <- function(api, pageSize = 1000, pageNumber = 0, ...) {
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    .invoke_bind(api, "getAllGenesUsingGET", pageSize = pageSize,
        pageNumber = pageNumber, ...)
}

#' @rdname cBioPortal
#'
#' @section API Metadata:
#' * queryGeneTable - Get a table for only the `genes` or `genePanelId` of
#'   interest. Gene inputs are identified with the `by` argument
#'
#' @examples
#' queryGeneTable(api = cbio, by = "entrezGeneId", genes = 7157)
#'
#' @export
queryGeneTable <- function(
        api,
        by = c("entrezGeneId", "hugoGeneSymbol"),
        genes = NA_character_,
        genePanelId = NA_character_
) {
    all.na <- function(x) all(is.na(x))

    if (all.na(genes) && all.na(genePanelId))
        stop("Provide either 'genes' or 'genePanelId'")

    by <- match.arg(by)
    geneIdType <- switch(
        by, entrezGeneId = "ENTREZ_GENE_ID", hugoGeneSymbol = 'HUGO_GENE_SYMBOL'
    )

    if (!all.na(genes))
        .invoke_bind(api, "fetchGenesUsingPOST",
            geneIdType = geneIdType, geneIds = as.character(genes))
    else
        getGenePanel(api, genePanelId = genePanelId)
}

#' @rdname cBioPortal
#'
#' @section Genes:
#' * getDataByGenes - Download data for a number of genes within
#'   `molecularProfileId` indicators, optionally a `sampleListId` can be
#'   provided.
#'
#' @examples
#' getDataByGenes(
#'     api = cbio,
#'     studyId = "acc_tcga",
#'     genes = 1,
#'     by = "entrezGeneId",
#'     molecularProfileIds = "acc_tcga_rna_seq_v2_mrna",
#'     sampleListId = "acc_tcga_rna_seq_v2_mrna"
#' )
#'
#' @export
getDataByGenes <-
    function(api, studyId = NA_character_, genes = NA_character_,
        genePanelId = NA_character_, by = c("entrezGeneId", "hugoGeneSymbol"),
        molecularProfileIds = NULL, sampleListId = NULL, sampleIds = NULL, ...)
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")
    if (!is.null(sampleListId))
        sampleIds <- samplesInSampleLists(api, sampleListId)[[1L]]
    else if (is.null(sampleIds))
        sampleIds <- allSamples(api, studyId)[["sampleId"]]

    by <- match.arg(by)

    feats <- queryGeneTable(api, by, genes, genePanelId)

    molData <- fetchData(
        api = api,
        studyId = studyId,
        molecularProfileIds = molecularProfileIds,
        entrezGeneIds = feats[["entrezGeneId"]],
        sampleIds = sampleIds
    )
    lapply(
        molData,
        function(x) suppressMessages({
            dplyr::left_join(x, feats)
        })
    )
}
