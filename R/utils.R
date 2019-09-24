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

.getMixedData <- function(x, name.field) {
    samplesAsCols <- .samplesAsCols(x)
    if (!any(samplesAsCols)) { return(.biocExtract(x)) }

    annote <- x[, !samplesAsCols]

    x <- data.matrix(x[, samplesAsCols])

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


.nonuniquesymbols <- function(vect) {
    if (is.null(vect))
        return(FALSE)
    as.logical(anyDuplicated(vect))
}

# vectorized version of finding name fields
.findNameFields <- function(x, names.field) {
    names.results <- Filter(length, lapply(names.field, function(nf)
        RTCGAToolbox:::.findCol(x, nf)))
    name.fields <- unlist(names.results, use.names = FALSE)
    if (!length(name.fields))
        name.fields <- NULL
    name.fields
}

# Get the column name of the name field that has unique identifiers at every
# row
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

.check_fun <- function(api, x, endpoint, colname, ...) {
    all(x %in% .invoke_bind(api, endpoint, ...)[[colname]])
}

.barg <- function(type) {
    switch(type,
        studyId = expression(list(keyword = element))[[1]],
        NULL
    )
}

.checkIdValidity <- function(api, element, ename = c("studyId", "genePanelId",
        "molecularProfileId", "sampleListId"), ...) {
    ename <- match.arg(ename)
    args <- .barg(ename)
    args <- eval(args)
    ord <- endpoint_map[endpoint_map[["what"]] == ename, , drop = TRUE]
    .check_fun(api = api, x = element,
        endpoint = ord[["how"]], colname = ename, args)
}

.generateIdConvert <- function(longid, shortid) {
    filler <- TCGAutils:::.uniqueDelim(longid)

    idlist <- strsplit(longid, filler)
    lens <- unique(lengths(idlist))
    if (length(lens) > 1)
        warning("Inconsistent sample codes (e.g., ", longid[1], ")")

    poss <- seq_len(lens)

    recomb <- function(x, range, filler)
        paste0(x[range], collapse = filler)

    validPOS <- vapply(poss,
        function(stoppage) {
            first <- vapply(idlist,
                function(x) recomb(x, seq_len(stoppage), filler),
                character(1L)
            )
            if (identical(unique(nchar(first)), unique(first)))
                return(FALSE)
            all(shortid %in% first)
        }, logical(1L)
    )

    if (sum(validPOS) != 1L)
        stop("Could not determine valid ID position cut-off")

    nrange <- seq_len(which(validPOS))
    args <- as.pairlist(alist(id =))
    body <- substitute({
        vapply(strsplit(id, g),
            function(x) paste0(x[z], collapse = g), character(1L)
        )
    }, list(g = filler, z = nrange))
    eval(call("function", args, body))
}
