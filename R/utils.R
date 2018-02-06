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

.nonuniquesymbols <- function(vect) {
    as.logical(anyDuplicated(vect))
}

.getMixedData <- function(x) {
    annotecols <- !startsWith(names(x), "TCGA")
    if (all(annotecols)) {
        return(.biocExtract(x))
    }
    annote <- x[, annotecols]

    hugodata <- RTCGAToolbox:::.hasHugoInfo(x)
    genecol <- grepl("^gene", names(annote), ignore.case = TRUE)

    x <- data.matrix(x[, !annotecols])
    if (hugodata) {
        hugoname <- RTCGAToolbox:::.findCol(annote, "Hugo_Symbol")
        geneSymbols <- annote[[hugoname]]
        if (!.nonuniquesymbols(geneSymbols))
            rownames(x) <- geneSymbols
    } else if (any(genecol)) {
        if (sum(genecol) == 1L)
        genenames <- annote[, genecol]
        if (!.nonuniquesymbols(genenames))
            rownames(x) <- genenames
    }
    x <- RTCGAToolbox:::.standardizeBC(x)
    SummarizedExperiment::SummarizedExperiment(SimpleList(x), rowData = annote)
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
