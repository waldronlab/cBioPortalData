.biocExtract <- function(object, names.field, colnames) {
    hasRanged <- RTCGAToolbox:::.hasRangeNames(object)
    if (hasRanged) {
        RTCGAToolbox:::.convertRangeBioc(object, names.field = names.field)
    } else {
        RTCGAToolbox:::.makeSummarizedExperimentFromDataFrame(
            object, names.field = names.field, colnames = colnames)
    }
}

.setBuild <- function(x, annotation) {
    build <- RTCGAToolbox:::.hasInfo(annotation, "ncbibuild")
    if (build) {
        buildno <- unique(RTCGAToolbox:::.getBuild(annotation, "ncbibuild"))
        isbuild <- TCGAutils::correctBuild(buildno, "NCBI")
        Seqinfo::genome(x) <- isbuild
    }
    x
}

.hasGRangesCols <- function(x) {
    !anyNA(TCGAutils::findGRangesCols(names(x))[c("seqnames", "start", "end")])
}

## drop rows with missing coordinates up front so that the ranges and the
## assay stay parallel; warn rather than drop silently
.dropNARanges <- function(x) {
    idx <- TCGAutils::findGRangesCols(names(x))[c("seqnames", "start", "end")]
    cols <- names(x)[idx]
    nas <- Reduce(`|`, lapply(cols, function(cl) is.na(x[[cl]])))
    if (any(nas)) {
        warning(
            "Dropped ", sum(nas), " of ", length(nas),
            " row(s) with missing genomic coordinates (",
            paste(cols, collapse = ", "), ")", call. = FALSE
        )
        x <- x[!nas, , drop = FALSE]
    }
    x
}

.makeGRanges <- function(x) {
    idx <- TCGAutils::findGRangesCols(names(x))
    fld <- function(f)
        if (is.na(idx[[f]])) NA_character_ else names(x)[idx[[f]]]
    strand <- fld("strand")
    noStrand <- is.na(strand) || all(is.na(x[[strand]]))
    gr <- GenomicRanges::makeGRangesFromDataFrame(
        x, keep.extra.columns = TRUE, seqnames.field = fld("seqnames"),
        start.field = fld("start"), end.field = fld("end"),
        strand.field = if (noStrand) "strand" else strand,
        ignore.strand = noStrand
    )
    ## row names come from the assay; leftover data.frame row names (e.g.,
    ## after dropping rows) would otherwise conflict with them
    names(gr) <- NULL
    gr
}

.getMutationData <- function(x, row.field) {
    x <- .dropNARanges(x)
    ridx <- na.omit(TCGAutils::findGRangesCols(names(x)))
    ranged <- x[, ridx]
    rowranges <- .makeGRanges(ranged)
    rowranges <- .setBuild(rowranges, x)

    others <- match(c("ncbibuild", "entrezgeneid", "hugogenesymbol"),
        tolower(names(x)))
    excl <- na.omit(c(ridx, others))

    rownames <- x[[row.field]]
    x <- as.matrix(x[, -excl])
    rownames(x) <- rownames

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

## GISTIC files require only peak_start / peak_end; region_start / region_end
## are optional extras and are preferred when present
.getGisticData <- function(x) {
    pfx <- if (all(c("region_start", "region_end") %in% names(x)))
        "region" else "peak"
    GenomicRanges::makeGRangesFromDataFrame(
        x, keep.extra.columns = TRUE, seqnames.field = "chromosome",
        start.field = paste0(pfx, "_start"), end.field = paste0(pfx, "_end"),
        ignore.strand = TRUE
    )
}

.getMixedData <- function(x, name.field) {
    samplesAsCols <- .TCGAsamplesAsCols(x)
    if (!length(x)) { return(x) }
    if (!any(samplesAsCols)) { return(.getcbiodata(x, name.field)) }

    annote <- x[, !samplesAsCols, drop = FALSE]
    hasRanged <- .hasGRangesCols(annote)
    if (hasRanged) {
        x <- .dropNARanges(x)
        annote <- x[, !samplesAsCols, drop = FALSE]
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

    ## build ranges only after all row filtering so they stay parallel to x
    if (hasRanged) {
        rowranges <- .makeGRanges(annote)
        rowranges <- .setBuild(rowranges, annote)
    }

    x <- RTCGAToolbox:::.standardizeBC(x)
    args <- list(assays = SimpleList(x))
    if (hasRanged)
        args <- append(args, list(rowRanges = rowranges))
    else
        args <- append(args, list(rowData = annote))
    do.call(SummarizedExperiment::SummarizedExperiment, args)
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
        if ("null" %in% strandvec) {
            x[strandvec %in% "null", strandname] <- "*"
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

.TCGAsamplesAsCols <- function(x) {
    startsWith(names(x), "TCGA")
}

.findColnames <- function(x, patients, columns) {
    dnames <- names(x)
    ptid_as_colnames <- dnames[dnames %in% patients]
    sampid_as_colnames <- dnames[dnames %in% columns]
    if (length(ptid_as_colnames))
        ptid_as_colnames
    else
        sampid_as_colnames
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
            message(msg)
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

.invoke_fun <- function(api, name, ...) {
    if (!is(api, "cBioPortal"))
        stop("<internal> Provide a 'cBioPortal' class API object")
    args <- Filter(Negate(is.null), list(...))
    api_call <- do.call(`$`, list(api, name))
    res <- do.call(api_call, args)
    msg <- httr::content(res)[["message"]]
    if (length(msg))
        message(msg)
    httr::stop_for_status(res)
    res
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_bind <- function(api, name, ...)  {
    .bind_content(.invoke_fun(api, name, ...))
}

