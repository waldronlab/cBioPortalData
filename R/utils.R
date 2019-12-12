utils::globalVariables("element")

.biocExtract <- function(object, names.field) {
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
                object, build = build, names.field = names.field)
    } else
        object <- RTCGAToolbox:::.makeSummarizedExperimentFromDataFrame(object,
            names.field = names.field)
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

    x <- as.matrix(x[, samplesAsCols])

    if (!is.null(name.field)) {
        rnames <- annote[[name.field]]
        nas <- is.na(rnames)
        ## remove NA values from both
        x <- x[!nas, ]
        annote <- annote[!nas, ]
        rownames(x) <- rnames[!nas]
        rownames(annote) <- rnames[!nas]
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
        x <- readr::type_convert(x)
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
        all(startsWith(col, "TCGA"))
    })
}

.hasMappers <- function(coldat) {
    pt <- RTCGAToolbox:::.findCol(coldat, "PATIENT_ID")
    samp <- RTCGAToolbox:::.findCol(coldat, "SAMPLE_ID")
    length(pt) && length(samp)
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

.check_fun <- function(api, x, endpoint, colname, use_cache, ...) {
    all(x %in% .invoke_bind(api, endpoint, use_cache, ...)[[colname]])
}

.barg <- function(type) {
    switch(type,
        studyId = expression(list(keyword = element))[[1]],
        NULL
    )
}

.checkIdValidity <- function(api, element, ename = c("studyId", "genePanelId",
        "molecularProfileId", "sampleListId"), use_cache = TRUE, ...) {
    if (all(is.na(element))) return(FALSE)
    ename <- match.arg(ename)
    args <- .barg(ename)
    args <- eval(args)
    ord <- endpoint_map[endpoint_map[["what"]] == ename, , drop = TRUE]
    .check_fun(api = api, x = element,
        endpoint = ord[["how"]], colname = ename, use_cache = use_cache, args)
}

.generateIdConvert <- function(longid, shortid) {
    filler <- TCGAutils:::.uniqueDelim(longid)
    if (!nchar(filler))
        stop("No clear delimiter in sample identifiers")

    idlist <- strsplit(longid, filler)
    lens <- unique(lengths(idlist))
    if (length(lens) > 1)
        warning("Inconsistent sample codes:\n",
            paste(Biobase::selectSome(longid), collapse = "\n"),
        call. = FALSE)

    poss <- seq_len(max(lens))

    posMAT <- lapply(setNames(poss, poss),
        function(endindx) {
            recomb <- vapply(idlist, function(splitid) {
                paste0(splitid[seq_len(endindx)], collapse = filler)
                }, character(1L)
            )
            recomb %in% shortid
        }
    )

    hitmat <- do.call(cbind, posMAT)
    ends <- apply(hitmat, 1L, sum)

    args <- as.pairlist(alist(id =))
    if (identical(length(unique(ends)), 1L))
        body <- substitute({
            vapply(strsplit(id, g),
                function(x) paste0(x[z], collapse = g), character(1L)
            )
        }, list(g = filler, z = seq_len(unique(ends))))
    else
        body <- substitute({
            mapply(function(x, y) paste0(x[seq_len(y)], collapse = g),
                strsplit(id, g), z)
        }, list(g = filler, z = ends))

    eval(call("function", args, body))
}

.invoke_fun <- function(api, name, use_cache = FALSE, ...) {
    if (!is(api, "cBioPortal"))
        stop("Provide a 'cBioPortal' class API object")
    ops <- names(AnVIL::operations(api))
    if (!name %in% ops)
        stop("<internal> operation name not found in API")

    if (use_cache) {
        .dollarCache(list(api, name), ...)
    } else {
        do.call(`$`, list(api, name))(...)
    }
}

.dollarCache <- function(appname, ...) {
    if (!is.list(appname))
        stop("<internal> Provide a list input as 'api$name'")
    digi <- digest::digest(list(appname, ...))
    loc <- .getHashCache(digi)
    if (file.exists(loc)) {
        load(loc)
    } else {
        op <- do.call(`$`, appname)(...)
        save(op, file = loc)
    }
    op
}

.bind_content <- function(x) {
    dplyr::bind_rows(
        httr::content(x)
    )
}

.invoke_bind <- function(api, name, use_cache = FALSE, ...) {
    .bind_content(.invoke_fun(api, name, use_cache, ...))
}
