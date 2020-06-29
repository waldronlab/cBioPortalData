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
#' for download. Please refer to the
#' \href{http://cbioportal.org/data_sets.jsp}{website} for the full list of
#' datasets. Users who would like to obtain a more fine grain approach to
#' downloading data from the cBioPortal API, refer to the `cBioPortalData`
#' function.
#'
#' @inheritParams downloadStudy
#'
#' @param split.field A character vector of possible column names for the column
#' that is used to identify samples in a mutations or copy number file.
#'
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
#'
#' @return A \linkS4class{MultiAssayExperiment} object
#'
#' @seealso \url{http://cbioportal.org/data_sets.jsp}, \link{cBioPortalData}
#'
#' @author Levi Waldron, M. Ramos
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
#' mae <- cBioDataPack("acc_tcga")
#'
#' @export
cBioDataPack <- function(cancer_study_id, use_cache = TRUE,
    split.field = c("Tumor_Sample_Barcode", "ID"),
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene")) {

    cancer_file <- downloadStudy(cancer_study_id, use_cache)

    exarg <- if (identical(.Platform$OS.type, "unix") &&
        Sys.info()["sysname"] != "Darwin")
        "--warning=no-unknown-keyword" else NULL

    filelist <- untar(cancer_file, list = TRUE, extras = exarg)
    filelist <- gsub("^\\.\\/", "", filelist)
    filekeepind <- grep("^\\._", basename(filelist), invert = TRUE)
    filelist <- filelist[filekeepind]
    ## Remove files that are corrupt / hidden (start with ._)
    datafiles <- grep(x = filelist, pattern = "data.*\\.(txt|seg)$",
        value = TRUE)
    datafiles <- c(datafiles, grep("meta_study", filelist, value = TRUE),
        grep("/LICENSE", filelist, value = TRUE))

    worktemp <- tempdir()
    untar(cancer_file, files = datafiles, exdir = worktemp,
        extras = exarg)

    exptfiles <- file.path(worktemp,
        grep("clinical|study|LICENSE|fusion|gistic", datafiles, invert = TRUE,
            value = TRUE))
    clinicalfiles <- file.path(worktemp,
        grep("clinical", datafiles, value = TRUE))
    mdatafile <- file.path(worktemp,
        grep("meta_study", datafiles, value = TRUE))
    licensefile <- file.path(worktemp,
        grep("/LICENSE", datafiles, value = TRUE))
    fusionExtra <- file.path(worktemp, grep("fusion", datafiles,
        value = TRUE, ignore.case = TRUE))
    gisticExtra <- file.path(worktemp, grep("gistic", datafiles,
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
            .hasMappers(readr::read_tsv(file, comment = "#", n_max = 5)),
            logical(1L)))
        if (length(clinwithcols) > 1) {
            clindatfile <- grep("sample|supp", names(clinwithcols),
                invert = TRUE, value = TRUE)
            if (length(clindatfile) > 1)
                clindatfile <- clindatfile[
                    which.max(vapply(clindatfile, function(file)
                    ncol(readr::read_tsv(file, n_max = 5L, comment = "#")),
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
    noTCGAstart <- !all(startsWith(idVector, "TCGA"))
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
