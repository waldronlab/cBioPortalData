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
        moldata <- getDataByGenePanel(api, by = by,
            genePanelId = genePanelId, studyId = studyId,
            molecularProfileId = molprof, sampleListId = sampleListId,
            check = check)
        if (grepl("mutation", molprof, ignore.case = TRUE))
            .getMutationData(moldata, by)
        else
            .getMixedData(moldata, by)
    })
    as(Filter(length, expers), "List")
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
        by = c("entrezGeneId", "hugoGeneSymbol"),
        idConvert = identity
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
    explist <- do.call(.portalExperiments, exargs)

    explist <- as(explist, "ExperimentList")
    clinargs <- match.args(clinicalData, call)
    clinargs <- eval.args(clinargs)
    clin <- do.call(clinicalData, clinargs)
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
