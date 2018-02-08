.hasRangeNames <- function(x) {
    if (is(x, "list")) { return(FALSE) }
    tcgacols <- all(grepl("^TCGA", names(x)))
    if (tcgacols) { return(FALSE) }
    if (!any(is.data.frame(x), is(x, "DataFrame"), is.matrix(x)))
        stop("(internal) 'x' must be rectangular")
    !all(is.na(TCGAutils::findGRangesCols(names(x))))
}

.biocExtract <- function(object, names.field) {
    hasRanged <- .hasRangeNames(object)
    build <- RTCGAToolbox:::.getBuild(object)
    if (hasRanged) {
        if (RTCGAToolbox:::.hasConsistentRanges(object))
            object <-
                RTCGAToolbox:::.makeRangedSummarizedExperimentFromDataFrame(
                    object, build = build, names.field = names.field)
        split.field <- RTCGAToolbox:::.findSampleCol(object)
        if (is.na(split.field) || !length(split.field))
            object <- RTCGAToolbox:::.makeGRangesFromDataFrame(object)
        else
            object <- RTCGAToolbox:::.makeRaggedExperimentFromDataFrame(
                object, build = build, names.field = names.field)
    } else
        object <- RTCGAToolbox:::.makeSummarizedExperimentFromDataFrame(object,
            names.field = names.field)
    return(object)
}


.getMixedData <- function(x, name.field) {
    annotecols <- !startsWith(names(x), "TCGA")
    if (all(annotecols)) {
        return(.biocExtract(x))
    }
    annote <- x[, annotecols]

    x <- data.matrix(x[, !annotecols])

    if (!is.null(name.field)) {
        rownames(x) <- annote[[name.field]]
    }

    x <- RTCGAToolbox:::.standardizeBC(x)
    SummarizedExperiment::SummarizedExperiment(SimpleList(x), rowData = annote)
}

.cleanHugo <- function(x) {
    hugodata <- RTCGAToolbox:::.hasInfo(x, "Hugo_Symbol")
    if (hugodata) {
        hugoname <- RTCGAToolbox:::.findCol(x, "Hugo_Symbol")
        compref <- grep("Composite.Element.REF", x[, hugoname, drop = TRUE],
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

.nonuniquesymbols <- function(vect) {
    if (is.null(vect))
        return(FALSE)
    as.logical(anyDuplicated(vect))
}

.findNameFields <- function(x, names.field) {
    names.results <- Filter(length, lapply(names.field, function(nf)
        RTCGAToolbox:::.findCol(x, nf)))
    name.fields <- unlist(names.results, use.names = FALSE)
    if (!length(name.fields))
        name.fields <- NULL
    name.fields
}

.getNameField <- function(x, names.field) {
    names.fields <- .findNameFields(x, names.field = names.field)
    vnames <- vapply(names.fields, function(ids) {
            !.nonuniquesymbols(x[[ids]])
        }, logical(1L))
    rname <- which(vnames)
    if (length(rname)) {
        return(names.fields[rname[[1L]]])
    }
    NULL
}
