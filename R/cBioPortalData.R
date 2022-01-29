.portalExperiments <- function(
    api, by, genePanelId, genes, studyId, molecularProfileIds,
    sampleListId, sampleIds
) {
    expers <- getDataByGenes(api, genes = genes, genePanelId = genePanelId,
        studyId = studyId, molecularProfileIds = molecularProfileIds,
        sampleListId = sampleListId, sampleIds = sampleIds, by = by)

    sampmap <- lapply(expers, function(x) {
        if (length(x)) {
            smap <- x[, c("molecularProfileId", "patientId", "sampleId")]
            names(smap) <- c("assay", "primary", "colname")
            smap
        } else {
            tibble::tibble(assay = character(0L), primary = character(0L),
                colname = character(0L))
        }
    })
    sampleMap <- dplyr::bind_rows(sampmap)

    experlist <- lapply(setNames(nm = names(expers)),
        function(molprof) {
            byGene <- expers[[molprof]]
            isMut <- grepl("mutation", molprof, ignore.case = TRUE)
            if (isMut)
                colsOI <- c(by, "chr", "startPosition", "endPosition",
                    "ncbiBuild", "sampleId", "mutationType")
            else
                colsOI <- c(by, "sampleId", "value")
            if (length(byGene)) {
                colsoi <- colsOI[colsOI %in% names(byGene)]

                if (isMut) {
                    res <- tidyr::pivot_wider(byGene[, colsoi],
                        names_from = "sampleId",
                        values_from = "mutationType",
                        values_fn = list(mutationType =
                            function(x) paste0(x, collapse = ";")
                        )
                    )
                    .getMutationData(res, by)
                } else {
                    res <- tidyr::pivot_wider(byGene[, colsoi],
                        names_from = "sampleId",
                        values_from = "value"
                    )
                    .getMixedData(res, by)
                }
            } else {
                SummarizedExperiment::SummarizedExperiment()
            }
        }
    )
    experlist <- as(Filter(length, experlist), "List")

    isTCGA <- grepl("tcga", studyId, ignore.case = TRUE)

    metalist <- lapply(names(experlist), function(molprof) {
        isMut <- grepl("mutation", molprof, ignore.case = TRUE)
        byGene <- expers[[molprof]]
        if (isMut) {
            colsOI <- c(by, "chr", "startPosition", "endPosition",
                "ncbiBuild", "sampleId", "mutationType")
            metaGene <- byGene[, !names(byGene) %in% colsOI]
            if (isTCGA) {
                ragex <- .getRagEx(byGene)
                experlist <<- c(experlist, list(proteinPos = ragex))
            }
        } else {
            colsOI <- c(by, "sampleId", "value")
            metaGene <- byGene[, !names(byGene) %in% colsOI]
            experlist <<- .updateRowData(metaGene, by, molprof, experlist)
        }
        metaGene
    })

    list(
        sampleMap = as(sampleMap, "DataFrame"),
        experiments = experlist,
        metadata = metalist
    )
}

.getRagEx <- function(metainfo, byname, assayname, explist) {
    ncbiBuildCol <- grep("ncbiBuild", names(metainfo), value = TRUE,
        fixed = TRUE)
    rangeSet <- cbind(
        metainfo[["proteinPosStart"]], metainfo[["proteinPosEnd"]]
    )
    metainfo[["proteinPosStart"]] <- apply(rangeSet, 1L, min, na.rm = TRUE)
    metainfo[["proteinPosEnd"]] <- apply(rangeSet, 1L, max, na.rm = TRUE)
    splitframe <- GenomicRanges::makeGRangesListFromDataFrame(metainfo,
        split.field = "sampleId", start.field = "proteinPosStart",
        end.field = "proteinPosEnd", keep.extra.columns = TRUE)
    ptIds <- TCGAutils::TCGAbarcode(names(splitframe))
    rex <- RaggedExperiment::RaggedExperiment(splitframe, colData =
        S4Vectors::DataFrame(row.names = ptIds)
    )
    if (length(ncbiBuildCol))
        genome(rex) <- TCGAutils::uniformBuilds(metainfo[[ncbiBuildCol]])
    rex
}

.updateRowData <- function(metainfo, byname, assayname, explist) {
    newby <- grep(byname, x = c("hugoGeneSymbol", "entrezGeneId"), value = TRUE,
        fixed = TRUE, invert = TRUE)
    stopifnot(is.character(newby), length(newby) == 1L)
    if (length(newby)) {
        exptoupdate <- explist[[assayname]]
        altNames <- unique(metainfo[[newby]])
        allName <- if (!is.null(altNames)) {
            all.equal(
                altNames,
                Reduce(intersect, split(metainfo[[newby]], metainfo[["patientId"]]))
            )
        } else { FALSE }
        if (isTRUE(allName)) {
            altDF <- DataFrame(altNames)
            names(altDF) <- newby
            rowData(exptoupdate) <- altDF
            explist[[assayname]] <- exptoupdate
        }
    }
    explist
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
    toeval <- !names(args) %in% c("api", "idConvert", "studyId")
    evalargs <- lapply(args[toeval], eval)
    stud <- dynGet("studyId")
    args["studyId"] <- stud
    args[toeval] <- evalargs
    args
}

update.args <- function(args) {
    molecularProfileIds <- args[["molecularProfileIds"]]
    api <- eval(args[["api"]])
    args[["api"]] <- api
    studyId <- args[["studyId"]]
    if (is.null(molecularProfileIds)) {
        molProfs <- molecularProfiles(api, studyId)
        ## data type not working yet
        ## https://github.com/cBioPortal/cbioportal/issues/7816
        args[["molecularProfileIds"]] <- molProfs[
            molProfs[["molecularAlterationType"]] != "STRUCTURAL_VARIANT",
            "molecularProfileId",
            drop = TRUE
        ]
    }

    args[["molecularProfileIds"]] <-
        setNames(nm = args[["molecularProfileIds"]])

    args
}

#' Download data from the cBioPortal API
#'
#' Obtain a `MultiAssayExperiment` object for a particular gene panel,
#' `studyId`, `molecularProfileIds`, and `sampleListIds` combination. Default
#' `molecularProfileIds` and `sampleListIds` are set to NULL for including all
#' data. This option is best for users who wish to obtain a section of the
#' study data that pertains to a specific molecular profile and gene panel
#' combination. For users looking to download the entire study data as provided
#' by the \url{https://cbioportal.org/datasets}, refer to `cBioDataPack`.
#'
#' @details We are able to succesfully represent 98 percent of the study
#'     identifiers as `MultiAssayExperiment` objects as obtained via
#'     `cBioPortalData` with the `IMPACT341` `genePanelId` as the example
#'     gene panel. Datasets that currently fail to import
#'     can be seen in the `getStudies(..., buildReport = TRUE)` dataset
#'     under the `"api_build"` column.
#'     Note that changes to the cBioPortal API may affect this rate at any
#'     time. If you encounter any issues, please open a GitHub issue at the
#'     \url{https://github.com/waldronlab/cBioPortalData/issues/} page with
#'     a fully reproducible example.
#'
#' @inheritParams cBioPortal
#'
#' @param check_build logical(1L) Whether to check the build status of the
#'     `studyId` using an internal dataset. This argument should be set to
#'     `FALSE` if using alternative `hostnames`, e.g.,
#'     'pedcbioportal.kidsfirstdrc.org'
#'
#' @md
#'
#' @examples
#'
#' cbio <- cBioPortal()
#'
#' samps <- samplesInSampleLists(cbio, "acc_tcga_rppa")[[1]]
#'
#' getGenePanelMolecular(
#'     cbio, molecularProfileIds = c("acc_tcga_rppa", "acc_tcga_linear_CNA"),
#'     samps
#' )
#'
#' acc_tcga <- cBioPortalData(
#'     cbio, by = "hugoGeneSymbol",
#'     studyId = "acc_tcga",
#'     genePanelId = "AmpliSeq",
#'     molecularProfileIds =
#'         c("acc_tcga_rppa", "acc_tcga_linear_CNA", "acc_tcga_mutations")
#' )
#'
#' @return A \linkS4class{MultiAssayExperiment} object
#'
#' @seealso \link{cBioDataPack}
#'
#' @export
cBioPortalData <-
    function(api, studyId = NA_character_,
        genePanelId = NA_character_,
        genes = NA_character_,
        molecularProfileIds = NULL,
        sampleListId = NULL,
        sampleIds = NULL,
        by = c("entrezGeneId", "hugoGeneSymbol"),
        check_build = TRUE
    )
{
    if (missing(api))
        stop("Provide a valid 'api' from 'cBioPortal()'")

    by <- match.arg(by)

    formals <- formals()
    formals[["by"]] <- by
    call <- std.args(match.call(), formals)
    exargs <- match.args(.portalExperiments, call)
    exargs <- eval.args(exargs)
    exargs <- update.args(exargs)

    if (check_build)
        .check_study_id_building(exargs[["studyId"]], "api_build")

    lists <- do.call(.portalExperiments, exargs)

    clinargs <- match.args(clinicalData, call)
    clinargs <- eval.args(clinargs)
    clin <- do.call(clinicalData, clinargs)
    clin <- as(clin, "DataFrame")

    # resolve duplicate IDs
    if (anyDuplicated(clin[["patientId"]])) {
        mets <- clin[duplicated(clin[["patientId"]]), ]
        metadata(clin) <- list(duplicated = mets)
        clin <- clin[!duplicated(clin[["patientId"]]), ]
    }
    rownames(clin) <- clin[["patientId"]]

    lists[["colData"]] <- clin
    do.call(MultiAssayExperiment, lists)
}
