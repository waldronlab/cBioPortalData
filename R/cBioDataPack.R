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
    }
    c(mdl, if (exists("lic")) LICENSE = lic)
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

cbioportal2clinicaldf <- function(file) {
    clin <- readr::read_tsv(file, comment = "#")
    clinmeta <- readr::read_tsv(file, col_names = FALSE, n_max = 2)
    clinmeta <- t(clinmeta)
    clinmeta <- sub("^\\#", "", clinmeta)
    colnames(clinmeta) <- c("column", "definition")
    clinmeta <- lapply(seq_along(colnames(clin)), function(i) {
        clinmeta[i, ]
    })
    names(clinmeta) <- colnames(clin)
    clin <- DataFrame(clin)
    metadata(clin) <- clinmeta
    clin <- .subBCLetters(clin)
    rownames(clin) <- clin[["PATIENT_ID"]]
    clin
}

.validStudyID <- function(cancer_study_id) {

    if (missing(cancer_study_id))
        stop("Provide a valid 'cancer_study_id' from 'studiesTable'")

    stopifnot(is.character(cancer_study_id),
        !is.na(cancer_study_id), length(cancer_study_id) == 1L)

    cancer_study_id <- tolower(cancer_study_id)
    ## Load dataset to envir
    loc_data <- new.env(parent = emptyenv())
    data("studiesTable", envir = loc_data, package = "cBioPortalData")
    studiesTable <- loc_data[["studiesTable"]]

    ## Ensure study ID is valid
    inTable <- cancer_study_id %in% studiesTable[["cancer_study_id"]]

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
#' from the `studiesTable` and retrieve the study tarball from cBioPortal.
#' These functions are used by `cBioDataPack` under the hood to download,
#' untar, and load the tarball datasets with caching. As stated in
#' `?cBioDataPack`, not all studies are currently working as
#' `MultiAssayExperiment` objects. As of July 2020, about ~80% of
#' datasets can be successfully imported into the `MultiAssayExperiment` data
#' class. Please open an issue if you would like the team to prioritize a
#' study. You may also check `studiesTable$pack_build` for a more current
#' status.
#'
#' @param cancer_study_id character(1) The study identifier from cBioPortal as
#' in \url{https://cbioportal.org/webAPI}
#'
#' @param use_cache logical(1) (default TRUE) create the default cache location
#' and use it to track downloaded data. If data found in the cache, data will
#' not be re-downloaded. A path can also be provided to data cache location.
#'
#' @param force logical(1) (default FALSE) whether to force re-download data from
#' remote location
#'
#' @param url_location character(1)
#' (default "https://cbioportal-datahub.s3.amazonaws.com") the URL location for
#' downloading packaged data. Can be set using the 'cBio_URL' option (see
#' `?cBioDataPack` for more details)
#'
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
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
    url_location = getOption("cBio_URL", .url_location))
{
    .validStudyID(cancer_study_id)

    url_file <- file.path(url_location, paste0(cancer_study_id, ".tar.gz"))

    if (is.character(use_cache) && length(use_cache) == 1L)
        cBioCache(directory = use_cache)
    else if (isTRUE(use_cache))
        cBioCache()
    else
        stop("Use 'setCache' or specify a download location")

    tryCatch({
        .download_data_file(url_file, cancer_study_id, verbose = TRUE,
            force = force)
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

#' @rdname downloadStudy
#'
#' @export
loadStudy <-
    function(
        filepath, names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene")
    )
{
    datafiles <- getRelevantFilesFromStudy(
        list.files(filepath, recursive = TRUE)
    )

    exptfiles <- file.path(filepath,
        grep("clinical|study|LICENSE|fusion|gistic", datafiles, invert = TRUE,
            value = TRUE))
    clinicalfiles <- file.path(filepath,
        grep("clinical", datafiles, value = TRUE))
    mdatafile <- file.path(filepath,
        grep("meta_study", datafiles, value = TRUE))
    licensefile <- file.path(filepath,
        grep("/LICENSE", datafiles, value = TRUE))
    fusionExtra <- file.path(filepath, grep("fusion", datafiles,
        value = TRUE, ignore.case = TRUE))
    gisticExtra <- file.path(filepath, grep("gistic", datafiles,
        value = TRUE, ignore.case = TRUE))

    expnames <- sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))
    expseq <- seq_along(exptfiles)
    names(expseq) <- expnames

    exptlist <- lapply(expseq, function(i, files, xpnames) {
        fname <- files[[i]]
        message(paste0("Working on: ", fname))
        dat <- as.data.frame(
            readr::read_tsv(fname, comment = "#"),
            check.names = FALSE)
        dat <- .cleanHugo(dat)
        dat <- .cleanStrands(dat)
        dat <- .standardizeBuilds(dat)

        names.field <- .findValidNames(dat, names.field)
        names.field <- .findUniqueField(dat, names.field)
        names.field <- .findMinDupField(dat, names.field)

        dat <- as(dat, "DataFrame")
        if (!RTCGAToolbox:::.hasExperimentData(dat))
            return(dat)
        cexp <- xpnames[[i]]
        if (grepl("meth", cexp)) {
            .getMixedData(dat, names.field)
        } else {
            .biocExtract(dat, names.field)
        }
    }, files = exptfiles, xpnames = expnames)

    names(exptlist) <-
        sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))

    .checkNonExpData <- function(exp) {
        is(exp, "GRanges") || is(exp, "DataFrame")
    }

    metadats <- Filter(.checkNonExpData, exptlist)
    exptlist <- Filter(function(expt) {!.checkNonExpData(expt)}, exptlist)

    if (length(clinicalfiles) > 1) {
        clinwithcols <- which(vapply(clinicalfiles, function(file)
            .hasMappers(
                readr::read_tsv(
                    file, comment = "#", n_max = 5, progress = FALSE
                )
            ),
            logical(1L)))
        if (length(clinwithcols) > 1) {
            clindatfile <- grep("sample|supp", names(clinwithcols),
                invert = TRUE, value = TRUE)
            if (length(clindatfile) > 1)
                clindatfile <- clindatfile[
                    which.max(vapply(clindatfile, function(file)
                    ncol(readr::read_tsv(
                        file, n_max = 5L, comment = "#", progress = FALSE
                    )),
                    integer(1L)))]
        } else
            clindatfile <- names(clinwithcols)
    } else {
        clindatfile <- clinicalfiles
    }

    coldata <- cbioportal2clinicaldf(clindatfile)
    mdat <- cbioportal2metadata(mdatafile, licensefile)

    if (length(fusionExtra))
        fudat <- readr::read_tsv(fusionExtra, comment = "#")
    else
        fudat <- list()

    if (length(gisticExtra))
        gist <- lapply(gisticExtra, function(x) {
            gfile <- readr::read_tsv(x, comment = "#")
            .getGisticData(gfile)
        })
    else
        gist <- list()

    mdat <- c(mdat, metadats, fudat, gist)
    exptlist <- MultiAssayExperiment::ExperimentList(exptlist)

    if (any(.TCGAcols(coldata))) {
        gmap <- TCGAutils::generateMap(exptlist, coldata,
            TCGAutils::TCGAbarcode)
    } else if (.hasMappers(coldata)) {
        gmap <- TCGAutils::generateMap(exptlist, coldata,
            sampleCol = "SAMPLE_ID", patientCol = "PATIENT_ID")
    } else {
        stop("Experiment data could not be mapped to colData")
    }

    MultiAssayExperiment(experiments = exptlist,
        colData = coldata, sampleMap = gmap, metadata = mdat)
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
#' @details The list of datasets can be found in the `studiesTable` dataset
#' by doing `data("studiesTable")`. Some datasets may not be available
#' for download and are not guaranteed to be represented as MultiAssayExperiment
#' data objects. After taking a random sample of 100
#' (using \code{set.seed(1234)}), we were able to succesfully represent about
#' 76 percent of the study identifiers as MultiAssayExperiment objects. Please
#' refer to the #' \href{http://cbioportal.org/data_sets.jsp}{website} for the
#' full list of available datasets. Users who would like to prioritize
#' particular datasets should open GitHub issues at the URL in the `DESCRIPTION`
#' file. For a more fine-grained approach to downloading data from the
#' cBioPortal API, refer to the `cBioPortalData` function.
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
#'
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
#'
#' @param ask A logical vector of length one indicating whether to prompt the
#' the user before downloading and loading study `MultiAssayExperiment`. If
#' TRUE, the user will be prompted to continue for studies that are not
#' currently building as `MultiAssayExperiment` based on previous testing
#' (in a non-interactive session, no data will be downloaded and built unless
#' `ask = FALSE`).
#'
#' @return A \linkS4class{MultiAssayExperiment} object
#'
#' @seealso \url{https://www.cbioportal.org/datasets}, \link{cBioPortalData}
#'
#' @author Levi Waldron, Marcel R., Ino dB.
#' @include utils.R
#'
#' @md
#'
#' @examples
#'
#' data(studiesTable)
#'
#' head(studiesTable[["cancer_study_id"]])
#'
#' # ask=FALSE for non-interactive use
#' mae <- cBioDataPack("acc_tcga", ask = FALSE)
#'
#' @export
cBioDataPack <- function(cancer_study_id, use_cache = TRUE,
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"), ask = TRUE) {

    denv <- new.env(parent = emptyenv())
    data("studiesTable", package = "cBioPortalData", envir = denv)
    studiesTable <- denv[["studiesTable"]]

    intable <- studiesTable[["cancer_study_id"]] %in% cancer_study_id
    if (!any(intable))
        stop("'cancer_study_id', ", cancer_study_id, ", not found.",
            " See 'data(\"studiesTable\")'.")

    builds <- studiesTable[["pack_build"]]
    hasbuilt <- unlist(builds[intable])

    if (!hasbuilt && any(builds)) {
        qtxt <- sprintf(
            paste0("Based on our tests, '%s' is not currently building.",
                "\n Proceed anyway? [y/n]: "),
            cancer_study_id
        )
        if (ask && .getAnswer(qtxt, allowed = c("y", "Y", "n", "N")) == "n")
            stop("'", cancer_study_id, "' is not yet supported.",
                " \n Use 'downloadStudy()' to obtain the study files.")
    }

    cancer_study_file <- downloadStudy(cancer_study_id, use_cache)
    exdir <- untarStudy(cancer_study_file)
    loadStudy(exdir, names.field)
}

