# previously http://download.cbioportal.org
.url_location <- "https://cbioportal-datahub.s3.amazonaws.com"

getRelevantFilesFromStudy <- function(filelist) {
    ## Remove files that are corrupt / hidden (start with ._)
    datafiles <- grep(x = filelist, pattern = "data.*\\.(txt|seg)$",
        value = TRUE)
    datafiles <- c(datafiles, grep("meta_study", filelist, value = TRUE),
        grep("/LICENSE", filelist, value = TRUE))
    datafiles
}

cbioportal2metadata <- function(meta_file, lic_file) {
    if (!length(meta_file) & !length(lic_file))
        return(list())
    md <- readLines(meta_file, warn = FALSE)
    mdl <- lapply(seq_along(md), function(i) {
        sub(".+: ", "", md[[i]])
    })
    names(mdl) <- sub(":.+", "", md)
    if (length(lic_file)) {
        lic <- readLines(lic_file, warn = FALSE)
        lic <- paste0(lic[lic != ""], collapse = "\n")
        lic <- list(LICENSE = lic)
    }
    c(mdl, if (exists("lic")) lic)
}

.subBCLetters <- function(df, ptID = "PATIENT_ID") {
    idVector <- df[[ptID]]
    allBC <- all(grepl("[A-Z]{4}.[0-9]{2}.[0-9]{4}", idVector))
    noTCGAstart <- is.character(idVector) && !all(startsWith(idVector, "TCGA"))
    if (allBC && noTCGAstart) {
        idVector <- gsub("^[A-Z]{4}", "TCGA", idVector)
        df[[ptID]] <- idVector
    }
    df
}

.silentRead <- function(file, comm = "#", mxlines = Inf, ...) {
    suppressMessages({
        readr::read_tsv(
            file, comment = comm, n_max = mxlines, progress = FALSE, ...
        )
    })
}

.processMeta <- function(clinmeta) {
    cnames <- unlist(unname(clinmeta[5L, ]))
    clinmeta <- clinmeta[-c(3L:5L), ]
    clinmeta <- t(clinmeta)
    clinmeta <- sub("^\\#", "", clinmeta)
    colnames(clinmeta) <- c("column", "definition")
    res <- lapply(setNames(seq_along(cnames), cnames), function(i) {
        clinmeta[i, ]
    })
    as(res, "DataFrame")
}

.getClinMeta <- function(clinfiles) {
    allmeta <- lapply(setNames(nm = clinfiles), function(x) {
        .silentRead(x, comm = "", mxlines = 5L, col_names = FALSE)
    })
    lapply(allmeta, .processMeta)
}

.readAll <- function(namedlist) {
    lapply(setNames(nm = names(namedlist)), function(x)
        .silentRead(x)
    )
}

.readSeparateMerge <- function(datalist) {
    alldata <- .readAll(datalist)
    Reduce(function(x, y) {
        merge(x, y, all = TRUE)
    }, alldata)
}

cbioportal2clinicaldf <- function(files) {
    if (length(files) > 1) {
        mappers <- lapply(setNames(nm = files), function(file)
            .whichMappers(.silentRead(file, mxlines = 5L))
        )
        hasMappers <- lengths(mappers) == 2L
        if (any(hasMappers)) {
            combdata <- mappers[hasMappers]
            clindata <- .readSeparateMerge(combdata)
        }
        ## try merge single mapper data to bigger merged
        singleCols <- lengths(mappers) == 1L
        if (all(singleCols)) {
            clindata <- .readSeparateMerge(mappers[singleCols])
        } else if (any(singleCols)) {
            singles <- .readAll(mappers[singleCols])
            clindata <- Reduce(function(x, y) {
                merge(x, y, all = TRUE)
            }, c(list(clindata), singles))
        }
    } else {
        clindata <- .silentRead(files, mxlines = 5L)
    }
    clinmeta <- .getClinMeta(files)
    clindata <- as(clindata, "DataFrame")
    metadata(clindata) <- clinmeta

    clindata <- .subBCLetters(clindata)
    rownames(clindata) <- clindata[["PATIENT_ID"]]
    clindata
}

.validStudyID <- function(cancer_study_id) {

    if (missing(cancer_study_id))
        stop("Provide a valid 'studyId' from 'getStudies'")

    stopifnot(is.character(cancer_study_id),
        !is.na(cancer_study_id), length(cancer_study_id) == 1L)

    cancer_study_id <- tolower(cancer_study_id)

    validStudies <- getStudies(cBioPortal())[["studyId"]]

    ## Ensure study ID is valid
    inTable <- cancer_study_id %in% validStudies

    if (!inTable)
        stop("Study identifier not found in look up table")
    else
        inTable
}

.download_data_file <-
    function(fileURL, cancer_study_id, verbose = FALSE, force = FALSE)
{
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (!length(rid)) {
        rid <- names(bfcadd(bfc, cancer_study_id, fileURL, download = FALSE))
    }
    if (!.cache_exists(bfc, cancer_study_id) || force) {
        if (verbose)
            message("Downloading study file: ", cancer_study_id, ".tar.gz")
            bfcdownload(bfc, rid, ask = FALSE)
    } else
        message("Study file in cache: ", cancer_study_id)

    bfcrpath(bfc, rids = rid)
}

.manageLocalFile <- function(cancer_study_id, inpath) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname", exact = TRUE)$rid
    if (!length(rid))
        stop("Can't update non-existing cache item")

    cachedir <- bfccache(bfc)
    finalname <- paste0(gsub("file", "", basename(tempfile())), "_",
        cancer_study_id, ".tar.gz")
    fileLoc <- file.path(cachedir, finalname)
    file.copy(inpath, fileLoc)

    bfcupdate(bfc, rids = rid, rpath = fileLoc)

    file.remove(inpath)

    bfcrpath(bfc, rids = rid)
}

.altDownload <- function(fileURL, cancer_study_id, verbose = FALSE) {
    if (verbose)
        message("Downloading study file: ", cancer_study_id, ".tar.gz")

    tmpFile <- file.path(tempdir(), paste0(cancer_study_id, ".tar.gz"))
    utils::download.file(fileURL, destfile = tmpFile, quiet = TRUE,
        method = "wget")

    .manageLocalFile(cancer_study_id, tmpFile)
}

#' @name downloadStudy
#'
#' @title Manually download, untar, and load study tarballs
#'
#' @description **Note** that these functions should be used when a particular
#' study is _not_ currently available as a `MultiAssayExperiment`
#' representation. Otherwise, use `cBioDataPack`. Provide a `cancer_study_id`
#' from `getStudies` and retrieve the study tarball from the cBio
#' Genomics Portal.  These functions are used by `cBioDataPack` under the hood
#' to download,untar, and load the tarball datasets with caching. As stated in
#' `?cBioDataPack`, not all studies are currently working as
#' `MultiAssayExperiment` objects. As of July 2020, about ~80% of
#' datasets can be successfully imported into the `MultiAssayExperiment` data
#' class. Please open an issue if you would like the team to prioritize a
#' study. You may also check `getStudies(buildReport = TRUE)$pack_build`
#' for the current status.
#'
#' @details When attempting to load a dataset using `loadStudy`, note that
#' the `cleanup` argument is set to `TRUE` by default. Change the argument
#' to `FALSE` if you would like to keep the untarred data in the `exdir`
#' location. `downloadStudy` and `untarStudy` are not affected by this change.
#' The tarball of the downloaded data is cached via `BiocFileCache` when
#' `use_cache` is `TRUE`.
#'
#' @param cancer_study_id character(1) The study identifier from cBioPortal as
#' in \url{https://cbioportal.org/webAPI}
#'
#' @param use_cache logical(1) (default TRUE) create the default cache location
#' and use it to track downloaded data. If data found in the cache, data will
#' not be re-downloaded. A path can also be provided to data cache location.
#'
#' @param ask logical(1) Whether to prompt the the user before downloading and
#'   loading study `MultiAssayExperiment` that is not currently building based
#'   on previous testing. Set to `interactive()` by default. In a
#'   non-interactive session, data download will be attempted; equivalent to
#'   `ask = FALSE`. The argument will also be used when a cache directory needs
#'   to be created when using `downloadStudy`.
#'
#' @param force logical(1) (default FALSE) whether to force re-download data
#' from remote location
#'
#' @param url_location character(1)
#' (default "https://cbioportal-datahub.s3.amazonaws.com") the URL location for
#' downloading packaged data. Can be set using the 'cBio_URL' option (see
#' `?cBioDataPack` for more details)
#'
#' @param names.field character() Possible column names for the
#' column that will used to label ranges for data such as mutations or copy
#' number (default: `c("Hugo_Symbol", "Entrez_Gene_Id", "Gene")`). Values are
#' cycled through and eliminated when no data present, or duplicates are found.
#' Values in the corresponding column must be unique in each row.
#'
#' @param cancer_study_file character(1) indicates the on-disk location
#' of the downloaded tarball
#'
#' @param exdir character(1) indicates the folder location to *put*
#' the contents of the tarball (default `tempdir()`; see also `?untar`)
#'
#' @param filepath character(1) indicates the folder location where
#' the contents of the tarball are *located* (usually the same as `exdir`)
#'
#' @param cleanup logical(1) whether to delete the `untar`-red contents from
#' the `exdir` folder (default TRUE)
#'
#' @return \itemize{
#'   \item {downloadStudy - The file location of the data tarball}
#'   \item {untarStudy - The directory location of the contents}
#'   \item {loadStudy - A MultiAssayExperiment-class object}
#' }
#'
#' @md
#'
#' @seealso \link{cBioDataPack}, \linkS4class{MultiAssayExperiment}
#'
#' @examples
#'
#' (acc_file <- downloadStudy("acc_tcga"))
#'
#' (file_dir <- untarStudy(acc_file, tempdir()))
#'
#' loadStudy(file_dir)
#'
#' @export
downloadStudy <- function(cancer_study_id, use_cache = TRUE, force = FALSE,
    url_location = getOption("cBio_URL", .url_location), ask = interactive())
{
    .validStudyID(cancer_study_id)

    url_file <- file.path(url_location, paste0(cancer_study_id, ".tar.gz"))

    if (is.character(use_cache) && length(use_cache) == 1L)
        cBioCache(directory = use_cache)
    else if (isTRUE(use_cache))
        cBioCache(ask = ask)
    else
        stop("Use 'setCache' or specify a download location")

    tryCatch(
        {
            .download_data_file(
                url_file, cancer_study_id, verbose = TRUE, force = force
            )
        },
        error = function(cond) {
            message("\n", cond)
            message("\nRetrying download with alternative function...")
            .altDownload(url_file, cancer_study_id, verbose = TRUE)
        }
    )
}

#' @rdname downloadStudy
#'
#' @export
untarStudy <- function(cancer_study_file, exdir = tempdir()) {
    exarg <- if (identical(.Platform$OS.type, "unix") &&
        Sys.info()["sysname"] != "Darwin")
        "--warning=no-unknown-keyword" else NULL

    filelist <- untar(cancer_study_file, list = TRUE, extras = exarg)
    filelist <- gsub("^\\.\\/", "", filelist)
    filekeepind <- grep("^\\._", basename(filelist), invert = TRUE)
    filelist <- filelist[filekeepind]
    datafiles <- getRelevantFilesFromStudy(filelist)

    folder <- basename(cancer_study_file)
    exdir <- file.path(exdir, gsub(".tar.gz", "", folder))
    if (!dir.exists(exdir))
        dir.create(exdir)

    untar(cancer_study_file, files = datafiles, exdir = exdir, extras = exarg)
    exdir
}

.preprocess_data <- function(file, exp_name, names.field, ptIDs) {
    if (is.null(exp_name))
        stop("<internal> 'exp_name' is NULL")

    message("Working on: ", file)
    dat <- utils::read.delim(
        file, sep = "\t", comment.char = "#", stringsAsFactors = FALSE,
        check.names = FALSE
    )
    dat <- .cleanHugo(dat)
    dat <- .cleanStrands(dat)
    dat <- .standardizeBuilds(dat)
    dat <- as(dat, "DataFrame")

    names.field <- .findValidNames(dat, names.field)
    names.field <- .findUniqueField(dat, names.field)
    names.field <- .findMinDupField(dat, names.field)

    tryCatch({
        if (!RTCGAToolbox:::.hasExperimentData(dat, ptIDs))
            dat
        else if (grepl("meth", exp_name, ignore.case = TRUE))
            .getMixedData(dat, names.field)
        else
            .biocExtract(dat, names.field, ptIDs)
    }, error = function(e) {
        err <- conditionMessage(e)
        warning(
            "Unable to import: ", exp_name, "\nReason: ", err, call. = FALSE
        )
        list()
    })
}

.loadExperimentsFromFiles <-
    function(fpath, dataFiles, names.field, colData)
{
    exptfiles <- file.path(fpath,
        grep("clinical|study|LICENSE|fusion|gistic", dataFiles, invert = TRUE,
            value = TRUE))
    expnames <- sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))
    names(exptfiles) <- expnames
    explist <- Map(
        function(x, y) {
            .preprocess_data(
                file = x, exp_name = y, names.field = names.field,
                ptIDs = colData[["PATIENT_ID"]]
            )
        },
        y = expnames, x = exptfiles
    )
    Filter(length, explist)
}

.isNonExpData <- function(exp) {
    is(exp, "GRanges") || is(exp, "DataFrame")
}

.readGISTIC <- function(filepath, datafiles, gist = list()) {
    gisticExtra <- .grepFiles("gistic", filepath, datafiles, ignore.case = TRUE)
    if (length(gisticExtra)) {
        gistics <- stats::setNames(gisticExtra, basename(gisticExtra))
        gist <- lapply(gistics, function(x) {
            gfile <- .silentRead(x)
            .getGisticData(gfile)
        })
    }
    gist
}

.readFUSION <- function(filepath, datafiles, fudat = list()) {
    fusionExtra <- .grepFiles("fusion", filepath, datafiles, ignore.case = TRUE)
    if (length(fusionExtra))
        fudat <- list(Fusion = .silentRead(fusionExtra))
    fudat
}

.grepFiles <- function(pattern, filepath, datafiles, ignore.case = FALSE) {
    file.path(
        filepath,
        grep(pattern, datafiles, value = TRUE, ignore.case = ignore.case)
    )
}

#' @rdname downloadStudy
#'
#' @export
loadStudy <- function(
    filepath,
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
    cleanup = TRUE
) {
    if (cleanup)
        on.exit(unlink(filepath, recursive = TRUE))

    datafiles <- getRelevantFilesFromStudy(
        list.files(filepath, recursive = TRUE)
    )

    mdatafile <- .grepFiles("meta_study", filepath, datafiles)
    licensefile <- .grepFiles("/LICENSE", filepath, datafiles)
    mdat <- cbioportal2metadata(mdatafile, licensefile)

    clinicalfiles <- .grepFiles("clinical", filepath, datafiles)
    coldata <- cbioportal2clinicaldf(clinicalfiles)

    explist <- .loadExperimentsFromFiles(
        fpath = filepath, dataFiles = datafiles,
        names.field = names.field, colData = coldata
    )

    slip <- split(explist, vapply(explist, .isNonExpData, logical(1L)))
    metadats <- slip[['TRUE']]
    explist <- MultiAssayExperiment::ExperimentList(slip[['FALSE']])

    fudat <- .readFUSION(filepath, datafiles)
    gist <- .readGISTIC(filepath, datafiles)

    mdat <- c(mdat, metadats, fudat, gist)

    if (any(.TCGAcols(coldata))) {
        gmap <- TCGAutils::generateMap(explist, coldata,
            TCGAutils::TCGAbarcode)
    } else if (.hasMappers(coldata)) {
        gmap <- TCGAutils::generateMap(explist, coldata,
            sampleCol = "SAMPLE_ID", patientCol = "PATIENT_ID")
    } else {
        stop("Experiment data could not be mapped to colData")
    }

    mdat <- c(mdat,
        unmapped = explist[names(explist) != unique(gmap[["assay"]])])

    MultiAssayExperiment(
        experiments = explist,
        colData = coldata,
        sampleMap = gmap,
        metadata = mdat
    )
}

.check_study_id_building <-
    function(
        cancer_study_id, build_type = c("pack_build", "api_build"), ask
    )
{
    build_type <- match.arg(build_type)
    denv <- .loadReportData()
    results <- denv[[build_type]]
    builds <- results[
        match(cancer_study_id, results[["studyId"]]), build_type
    ]

    if (is.na(builds))
        stop("'studyId', ", cancer_study_id, ", not found.",
            " See 'getStudies()'.")

    if (!builds) {
        qtxt <- sprintf(
            paste0(
                "'getStudies' reports that '%s' is not currently building.\n",
                "  Use 'downloadStudy()' to obtain the study data.\n",
                "  Proceed anyway? [y/n]: "
            ),
            cancer_study_id
        )
        if (ask && .getAnswer(qtxt, allowed = c("y", "Y", "n", "N")) == "n")
            stop("'", cancer_study_id, "' is not yet supported.")
    }

}

#' @name cBioDataPack
#'
#' @title Obtain pre-packaged data from cBioPortal and represent as
#' a MultiAssayExperiment object
#'
#' @description The `cBioDataPack` function allows the user to
#' download and process cancer study datasets found in MSKCC's cBioPortal.
#' Output datasets use the \linkS4class{MultiAssayExperiment} data
#' representation to faciliate analysis and data management operations.
#'
#' @details The full list of study identifiers (`studyId`s) can obtained from
#' `getStudies()`. Currently, only ~ 72% of datasets can be represented as
#' `MultiAssayExperiment` data objects from the data tarballs. Refer to
#' `getStudies(..., buildReport = TRUE)` and its `"pack_build"` column to see
#' which study identifiers are not building. Users who would like to prioritize
#' particular datasets should open GitHub issues at the URL in the
#' `DESCRIPTION` file. For a more fine-grained approach to downloading data
#' from the cBioPortal API, refer to the `cBioPortalData` function.
#'
#' @section cBio_URL:
#' The `cBioDataPack` function accesses data from the `cBio_URL` option.
#' By default, it points to an Amazon S3 bucket location. Previously, it
#' pointed to 'http://download.cbioportal.org'. This recent change
#' (> 2.1.17) should provide faster and more reliable downloads for all users.
#' See the URL using `cBioPortalData:::.url_location`. This can be changed
#' if there are mirrors that host this data by setting the `cBio_URL` option
#' with `getOption("cBio_URL", "https://some.url.com/")` before running the
#' function.
#'
#' @inheritParams downloadStudy
#' @inheritParams cBioPortalData
#'
#' @return A \linkS4class{MultiAssayExperiment} object
#'
#' @seealso \url{https://www.cbioportal.org/datasets}, \link{cBioPortalData},
#'   \link{removePackCache}
#'
#' @author Levi Waldron, Marcel R., Ino dB.
#' @include utils.R
#'
#' @md
#'
#' @examples
#'
#' cbio <- cBioPortal()
#'
#' head(getStudies(cbio)[["studyId"]])
#'
#' mae <- cBioDataPack("acc_tcga")
#'
#' @export
cBioDataPack <- function(cancer_study_id, use_cache = TRUE,
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"),
    cleanup = TRUE, ask = interactive(), check_build = TRUE)
{
    if (check_build)
        .check_study_id_building(cancer_study_id, "pack_build", ask = ask)

    cancer_study_file <- downloadStudy(
        cancer_study_id, use_cache = use_cache, ask = ask
    )
    exdir <- untarStudy(cancer_study_file)
    loadStudy(exdir, names.field, cleanup)
}

