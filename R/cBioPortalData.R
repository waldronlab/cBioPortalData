.portalExperiments <-
    function(api, by, genePanelId, studyId, molecularProfileIds, sampleListId,
        check)
{
    if (is.null(molecularProfileIds)) {
        molecularProfileIds <-
            molecularProfiles(api, studyId)[["molecularProfileId"]]
    } else { check <- TRUE }

    molecularProfileIds <- stats::setNames(molecularProfileIds,
        molecularProfileIds)

    expers <- lapply(molecularProfileIds, function(molprof) {
        getDataByGenePanel(api, by = by,
            genePanelId = genePanelId, studyId = studyId,
            molecularProfileId = molprof, sampleListId = sampleListId,
            check = check)
    })
    
    sampmap <- lapply(expers, function(x) {
        smap <- x[, c("molecularProfileId", "patientId", "sampleId")]
        names(smap) <- c("assay", "primary", "colname")
        smap
    })
    sampleMap <- dplyr::bind_rows(sampmap)
    
    explist <- lapply(molecularProfileIds, function(molprof) {
        isMut <- grepl("mutation", molprof, ignore.case = TRUE)
        byGene <- expers[[molprof]]
        if (isMut)
            colsOI <- c("entrezGeneId","chr", "startPosition", "endPosition",
                "ncbiBuild", "sampleId", "mutationType")
        else
            colsOI <- c("entrezGeneId", "sampleId", "value")
        colsoi <- colsOI[colsOI %in% names(byGene)]
        if (isMut) {
            res <- tidyr::pivot_wider(byGene[, colsoi], names_from = "sampleId",
                values_from = "mutationType",
                values_fn = list(mutationType =
                    function(x) paste0(x, collapse = ";")))
            .getMutationData(res, by)
        } else {
            res <- tidyr::pivot_wider(byGene[, colsoi], names_from = "sampleId",
                values_from = "value")
            .getMixedData(res, by)
        }
    })
    as(Filter(length, explist), "List")
    
    metalist <- lapply(names(expers), function(molprof) {
        isMut <- grepl("mutation", molprof, ignore.case = TRUE)
        byGene <- expers[[molprof]]
        if (isMut)
            colsOI <- c("entrezGeneId","chr", "startPosition", "endPosition",
                "ncbiBuild", "sampleId", "mutationType")
        else
            colsOI <- c("entrezGeneId", "sampleId", "value")
        byGene[, !names(byGene) %in% colsOI]
    })
    
    list(
        sampleMap = as(sampleMap, "DataFrame"),
        experiments = explist,
        metadata = metalist
    )
}

std.args <- function(call, formals) {
    callargs <- as.list(call)[-1]
    toadd <- setdiff(names(formals), names(callargs))
    call[toadd] <- formals[toadd]
    call
}

match.args <- function(fun, call, ...) {
    funfor <- formals(fun)
    exargs <- intersect(names(funfor), names(call))
    c(as.list(call)[-1][exargs], ...)
}

eval.args <- function(args) {
    toeval <- !names(args) %in% c("api", "idConvert")
    evalargs <- lapply(args[toeval], eval)
    args[toeval] <- evalargs
    args
}

.buildMap <- function(api, studyid, elist) {
    samptable <- allSamples(api, studyid)[, c("patientId", "sampleId")]
    cnames <- colnames(elist)
    smap <- lapply(cnames, function(cnms) {
        samptable[samptable[["sampleId"]] %in% cnms, ]
    })
    nmap <- cbind(rep(names(cnames), vapply(smap, nrow, integer(1L))),
        dplyr::bind_rows(smap))
    names(nmap) <- c("assay", "primary", "colname")
    nmap
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
#' cbio <- cBioPortal()
#'
#' cBioPortalData(cbio, by = "hugoGeneSymbol", studyId = "acc_tcga",
#'     genePanelId = "IMPACT341",
#'     molecularProfileIds = c("acc_tcga_rppa", "acc_tcga_linear_CNA")
#' )
#'
#' @export
cBioPortalData <-
    function(api, studyId = NA_character_,
        genePanelId = NA_character_,
        molecularProfileIds = NULL,
        sampleListId = NULL,
        by = c("entrezGeneId", "hugoGeneSymbol")
    )
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    validStudy <- .checkIdValidity(api, element = studyId, ename = "studyId")
    if (!validStudy)
        stop("Provide a valid 'studyId' from 'getStudies()'")

    validGP <-
        .checkIdValidity(api, element = genePanelId, ename = "genePanelId")
    if (!validGP)
        stop("Provide a valid 'genePanelId' from 'genePanels()'")

    by <- match.arg(by)

    formals <- formals()
    formals[["by"]] <- by
    call <- std.args(match.call(), formals)
    exargs <- match.args(.portalExperiments, call, check = FALSE)
    exargs <- eval.args(exargs)
    lists <- do.call(.portalExperiments, exargs)

    clinargs <- match.args(clinicalData, call)
    clinargs <- eval.args(clinargs)
    clin <- do.call(clinicalData, clinargs)
    clin <- as.data.frame(clin)
    rownames(clin) <- clin[["patientId"]]
    
    lists[["colData"]] <- clin
    
    if (isEmpty(lists[["sampleMap"]]))
        sampmap <- .buildMap(api, studyId, lists[["ExperimentList"]])
    #    sampmap <- TCGAutils::generateMap(experiments = explist,
    #        colData = clin, idConverter = idConvert)

    do.call(MultiAssayExperiment, lists)
}
