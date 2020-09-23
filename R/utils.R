utils::globalVariables("element")

.biocExtract <- function(object, names.field, colnames) {
    hasRanged <- RTCGAToolbox:::.hasRangeNames(object)
    build <- RTCGAToolbox:::.getBuild(object)
    if (hasRanged) {
        if (RTCGAToolbox:::.hasConsistentRanges(object))
            object <-
                RTCGAToolbox:::.makeRangedSummarizedExperimentFromDataFrame(
                    object, build = build, names.field = names.field)
        split.field <- RTCGAToolbox:::.findSampleCol(object)
        if (is.na(split.field) || !length(split.field))
            object <- RTCGAToolbox:::.makeGRangesFromDataFrame(object,
                build = build)
        else
            object <- RTCGAToolbox:::.makeRaggedExperimentFromDataFrame(
                object, build = build, names.field = names.field,
                split.field = split.field
            )
    } else {
        object <- RTCGAToolbox:::.makeSummarizedExperimentFromDataFrame(object,
            names.field = names.field, colnames = colnames)
    }
    return(object)
}

.getMutationData <- function(x, row.field) {
    build <- RTCGAToolbox:::.hasInfo(x, "ncbibuild")
    if (build)
        buildno <- RTCGAToolbox:::.getBuild(x)
    rownames <- x[[row.field]]
    ridx <- na.omit(TCGAutils::findGRangesCols(names(x)))
    ranged <- x[, ridx]
    others <- match(c("ncbibuild", "entrezgeneid", "hugogenesymbol"),
        tolower(names(x)))
    excl <- na.omit(c(ridx, others))
    x <- as.matrix(x[, -excl])
    rownames(x) <- rownames
    rowranges <- RTCGAToolbox:::.makeGRangesFromDataFrame(ranged)
    if (build)
        genome(rowranges) <- buildno
    SummarizedExperiment::SummarizedExperiment(assays = x,
        rowRanges = rowranges)
}

.getcbiodata <- function(x, name.field) {
    if (!is.null(name.field) || length(name.field)) {
        rowscol <- match(name.field, names(x))
        rnames <- x[[name.field]]
        x <- as.matrix(x[, -rowscol])
        rownames(x) <- rnames
    }
    SummarizedExperiment::SummarizedExperiment(x)
}

.getGisticData <- function(x) {
    RTCGAToolbox:::.makeGRangesFromDataFrame(x)
}

.getMixedData <- function(x, name.field) {
    samplesAsCols <- .samplesAsCols(x)
    if (!length(x)) { return(x) }
    if (!any(samplesAsCols)) { return(.getcbiodata(x, name.field)) }

    annote <- x[, !samplesAsCols, drop = FALSE]
    hasRanged <- RTCGAToolbox:::.hasRangeNames(annote)
    if (hasRanged) {
        rowranges <- RTCGAToolbox:::.makeGRangesFromDataFrame(annote)
        build <- RTCGAToolbox:::.hasInfo(annote, "ncbibuild")
        if (build) {
            isbuild <- RTCGAToolbox:::.getBuild(annote, "ncbibuild")

            genome(rowranges) <- isbuild
        }
    }

    x <- data.matrix(x[, samplesAsCols])

    if (!is.null(name.field)) {
        rnames <- annote[[name.field]]
        nasdu <- is.na(rnames) | duplicated(rnames)
        ## remove NA and duplicate values from both
        x <- x[!nasdu, , drop = FALSE]
        annote <- as.data.frame(annote[!nasdu, ])
        rownames(x) <- rnames[!nasdu]
        rownames(annote) <- rnames[!nasdu]
    }

    x <- RTCGAToolbox:::.standardizeBC(x)
    SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(x),
        if (hasRanged) rowRanges = rowranges else rowData = annote
    )
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
        x <- suppressMessages(readr::type_convert(x))
    }
    return(x)
}

.cleanStrands <- function(x) {
    strandData <- RTCGAToolbox:::.hasInfo(x, "Strand")
    if (strandData) {
        strandname <- RTCGAToolbox:::.findCol(x, "Strand")
        strandvec <- x[, strandname, drop = TRUE]
        if (any(c(1L, -1L) %in% strandvec)) {
            newStrand <- rep("*", length(strandvec))
            newStrand[strandvec %in% 1L] <- "+"
            newStrand[strandvec %in% -1L] <- "-"
            x[, strandname] <- newStrand
        }
    }
    return(x)
}

.standardizeBuilds <- function(x) {
    ncbi <- RTCGAToolbox:::.findCol(x, "NCBI_Build")
    if (length(ncbi)) {
    x[[ncbi]] <- TCGAutils::uniformBuilds(x[[ncbi]])
    }
    return(x)
}

# vectorized version of finding name fields
.findValidNames <- function(x, names.field) {
    names.results <- Filter(length, lapply(names.field, function(nf)
        RTCGAToolbox:::.findCol(x, nf)))
    name.fields <- unlist(names.results, use.names = FALSE)
    if (!length(name.fields))
        name.fields <- NULL
    name.fields
}

# Get the column name of the name field that has unique identifiers at every
# row
.findUniqueField <- function(x, names.fields) {
    if (length(names.fields) > 1L) {
        vnames <- vapply(names.fields, function(id) {
            if (is.null(x[[id]]))
                NULL
            else
                as.logical(anyDuplicated(x[[id]]))
            }, logical(1L))
        if (!length(vnames))
            NULL
        else
            names(vnames)[which(vnames)]
    } else {
        names.fields
    }
}

.findMinDupField <- function(x, names.fields) {
    if (length(names.fields) > 1L) {
        dupsum <- vapply(names.fields, function(id) {
            if (is.null(x[[id]]))
                Inf
            else
                sum(duplicated(x[[id]]))
        }, numeric(1L))
        vmin <- min(dupsum)
        if (all(is.infinite(vmin)))
            NULL
        else # take first when same number of duplicates
            names.fields[dupsum == vmin][[1L]]
    } else {
        names.fields
    }
}

.samplesAsCols <- function(x) {
    startsWith(names(x), "TCGA")
}

.TCGAcols <- function(df) {
    apply(df, 2L, function(col) {
        all(grepl("^TCGA", col, ignore.case = TRUE))
    })
}

.whichMappers <- function(coldat) {
    c(
        RTCGAToolbox:::.findCol(coldat, "PATIENT_ID"),
        RTCGAToolbox:::.findCol(coldat, "SAMPLE_ID")
    )
}

.hasMappers <- function(coldat) {
    mapps <- .whichMappers(coldat)
    length(mapps) == 2L
}

.getAnswer <- function(msg, allowed)
{
    if (interactive()) {
        repeat {
            cat(msg)
            answer <- readLines(n = 1)
            if (answer %in% allowed)
                break
        }
        tolower(answer)
    } else {
        "n"
    }
}

endpoint_map <- data.frame(
    what = c("studyId", "genePanelId", "molecularProfileId", "sampleListId"),
    how = c("getAllStudiesUsingGET", "getAllGenePanelsUsingGET",
        "getAllMolecularProfilesUsingGET", "getAllSampleListsUsingGET"),
    stringsAsFactors = FALSE
)

.dollarCache <- function(appname, ...) {
    if (!is.list(appname))
        stop("<internal> Provide a list input as 'api$name'")
    digi <- digest::digest(list(appname, ...))
    loc <- .getHashCache(digi)
    if (file.exists(loc)) {
        load(loc)
    } else {
        op <- do.call(`$`, appname)(...)
        save(op, file = loc, compress = FALSE)
    }
    op
}

.invoke_fun <- function(api, name, use_cache = FALSE, ...) {
    if (!is(api, "cBioPortal"))
        stop("Provide a 'cBioPortal' class API object")

    if (use_cache) {
        .dollarCache(list(api, name), ...)
    } else {
        do.call(`$`, list(api, name))(...)
    }
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_bind <- function(api, name, use_cache, ...)  {
    .bind_content(.invoke_fun(api, name, use_cache, ...))
}

