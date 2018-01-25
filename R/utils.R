.hasRangeNames <- function(x) {
    if (is(x, "list")) { return(FALSE) }
    tcgacols <- all(grepl("^TCGA", names(x)))
    if (tcgacols) { return(FALSE) }
    if (!any(is.data.frame(x), is(x, "DataFrame"), is.matrix(x)))
        stop("(internal) 'x' must be rectangular")
    !all(is.na(TCGAutils::findGRangesCols(names(x))))
}

.biocExtract <- function(object) {
    hasRanged <- .hasRangeNames(object)
    if (hasRanged) {
        if (RTCGAToolbox:::.hasBuildInfo(object)) {
            GBuild <- RTCGAToolbox:::.getBuild(object)
        }
        bld <- if (exists("GBuild")) { GBuild } else { NULL }
        if (RTCGAToolbox:::.hasConsistentRanges(object)) {
            object <-
            RTCGAToolbox:::.makeRangedSummarizedExperimentFromDataFrame(
                object, build = bld)
        }
        split.field <- RTCGAToolbox:::.findSampleCol(object)
        if (is.na(split.field) || !length(split.field)) {
            object <- RTCGAToolbox:::.makeGRangesFromDataFrame(object)
        } else {
            object <- RTCGAToolbox:::.makeRaggedExperimentFromDataFrame(
                object, build = bld)
        }
    } else {
        object <- RTCGAToolbox:::.standardizeBC(object)
        metadat <- metadata(object)
        object <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    return(object)
}

.getGISTIC <- function(x) {
    annoteCols <- !startsWith(names(x), "TCGA")
    if (all(annoteCols)) {
        return(.biocExtract(x))
    }
    annoteRowDF <- x[, annoteCols]
    genecol <- grepl("^gene", names(annoteRowDF), ignore.case = TRUE)
    rownames(annoteRowDF) <- annoteRowDF[, genecol]
    x <- x[, !annoteCols]
    x <- vapply(x, type.convert, numeric(nrow(x)))
    x <- RTCGAToolbox:::.standardizeBC(x)
    SummarizedExperiment::SummarizedExperiment(SimpleList(x), rowData = annoteRowDF)
}

.getMethyl <- function(x) {
    x <- as.data.frame(x)
    annote <- x[, !startsWith(names(x), "TCGA")]
    isNumRow <- all(grepl("^[0-9]*$",
        sample(rownames(x), size = 100L, replace = TRUE)))
    if (isNumRow) {
        geneSymbols <- annote[, grep("symbol", names(annote),
            ignore.case = TRUE, value = TRUE)]
        rNames <- geneSymbols
    } else { rNames <- rownames(x) }
    dm <- data.matrix(x[, startsWith(names(x), "TCGA")])
    rownames(dm) <- rNames
    dm <- RTCGAToolbox:::.standardizeBC(dm)
    SummarizedExperiment::SummarizedExperiment(SimpleList(dm), rowData = annote)
}
