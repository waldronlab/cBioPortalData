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
        if (RTCGAToolbox:::.hasBuildInfo(object))
            GBuild <- RTCGAToolbox:::.getBuild(object)
        bld <- if (exists("GBuild")) { GBuild } else { NULL }
        if (RTCGAToolbox:::.hasConsistentRanges(object))
            object <-
                RTCGAToolbox:::.makeRangedSummarizedExperimentFromDataFrame(
                    object, build = bld)
        split.field <- RTCGAToolbox:::.findSampleCol(object)
        if (is.na(split.field) || !length(split.field))
            object <- RTCGAToolbox:::.makeGRangesFromDataFrame(object)
        else
            object <- RTCGAToolbox:::.makeRaggedExperimentFromDataFrame(
                object, build = bld)
    } else
        object <- RTCGAToolbox:::.makeSummarizedExperimentFromDataFrame(object)
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

.cleanHugo <- function(x) {
    hugodata <- RTCGAToolbox:::.hasHugoInfo(x)
    if (hugodata) {
        hugoname <- RTCGAToolbox:::.findCol(x, "Hugo_Symbol")
        compref <- grep("Composite Element REF", x[, hugoname, drop = TRUE],
            ignore.case = TRUE)
        if (length(compref))
            x <- x[-compref, , drop = FALSE]
        hugos <- x[, hugoname, drop = TRUE]
        hugos <- vapply(strsplit(hugos, "|", TRUE), `[`, character(1L), 1L)
        x[, hugoname] <- hugos
        x <- readr::type_convert(x)
    }
    return(x)
}
