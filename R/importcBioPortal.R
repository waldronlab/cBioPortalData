download_data_file <- function(fileURL, cancer_study_id, verbose = FALSE) {
    bfc <- .get_cache()
    rid <- bfcquery(bfc, cancer_study_id, "rname")$rid
    if (!length(rid)) {
        if( verbose )
            message( "Downloading study file: ", cancer_study_id, ".tar.gz")
        rid <- names(bfcadd(bfc, cancer_study_id, fileURL))
    }
    if (!.cache_exists(bfc, cancer_study_id))
        bfcdownload(bfc, rid, ask = FALSE)

    bfcrpath(bfc, rids = rid)
}

#' @title Convert a data file downloaded from MSKCC's cBioPortal to
#' a MultiAssayExperiment object
#'
#'
#' @param cancer_study_id The cBioPortal study identifier
#' @param use_cache logical (default TRUE) create the default cache location
#' and use it to track downloaded data. If data found in the cache, data will
#' not be re-downloaded
#' @param split.field A character vector of possible column names for the column
#' that is used to identify samples in a mutations or copy number file.
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
#'
#' @return A \code{MultiAssayExperiment} object
#' @seealso \url{http://cbioportal.org/data_sets.jsp}
#'
#' @author Levi Waldron, M. Ramos
#' @include utils.R
#'
#' @examples
#'
#' data(studiesTable)
#'
#' (laml <- studiesTable[["cancer_study_id"]][3L])
#'
#' mae <- importcBioPortal(laml)
#'
#' @export importcBioPortal
importcBioPortal <- function(cancer_study_id, use_cache = TRUE,
    split.field = c("Tumor_Sample_Barcode", "ID"),
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id", "Gene")) {

    if (missing(cancer_study_id))
        stop("Provide a valid 'cancer_study_id' from 'studiesTable'")
    else
        cancer_study_id <- tolower(cancer_study_id)
    ## Load dataset to envir
    loc_data <- new.env(parent = emptyenv())
    data("studiesTable", envir = loc_data)
    studiesTable <- loc_data[["studiesTable"]]

    ## Ensure study ID is valid
    inTable <- cancer_study_id %in% studiesTable[["cancer_study_id"]]
    if (!S4Vectors::isSingleString(cancer_study_id) || !inTable)
        stop("Provide a single and valid study identifier")

    url_location <- paste0("https://media.githubusercontent.com/media/",
        "cBioPortal/datahub/master/public")
    url_file <- file.path(url_location, paste0(cancer_study_id, ".tar.gz"))

    if (use_cache)
        setCache(verbose = FALSE, ask = FALSE)
    else
        stop("Use 'setCache' to specify a download location")

    cancer_file <- download_data_file(url_file, cancer_study_id, verbose = TRUE)

    fileList <- untar(cancer_file, list = TRUE)
    ## Remove files that are corrupt / hidden (start with ._)
    datafiles <- grep(fileList, pattern = "[^\\._]data.+\\.(txt|seg)$", value = TRUE)
    datafiles <- c(datafiles, grep("[^\\._]meta_study", fileList, value = TRUE),
        grep("/LICENSE", fileList, value = TRUE))

    worktemp <- tempdir()
    untar(cancer_file, files = datafiles, exdir = worktemp)

    exptfiles <- file.path(worktemp,
        grep("clinical|study|LICENSE|fusion", datafiles, invert = TRUE,
            value = TRUE))
    clinicalfiles <- file.path(worktemp,
        grep("clinical", datafiles, value = TRUE))
    mdatafile <- file.path(worktemp,
        grep("meta_study", datafiles, value = TRUE))
    licensefile <- file.path(worktemp,
        grep("/LICENSE", datafiles, value = TRUE))
    fusionExtra <- file.path(worktemp, grep("fusion", datafiles,
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

        name.field <- .getNameField(dat, names.field = names.field)
        dat <- as(dat, "DataFrame")
        if (!RTCGAToolbox:::.hasExperimentData(dat))
            return(dat)
        cexp <- xpnames[[i]]
        if (grepl("meth", cexp) || grepl("gist", cexp)) {
            .getMixedData(dat, name.field)
        } else {
            .biocExtract(dat, name.field)
        }
    }, files = exptfiles, xpnames = expnames)

    names(exptlist) <-
        sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))

    .checkNonExpData <- function(exp) {
        is(exp, "GRanges") || is(exp, "DataFrame")
    }

    metadats <- Filter(.checkNonExpData, exptlist)
    exptlist <- Filter(function(expt) {!.checkNonExpData(expt)}, exptlist)

    clindatfile <- grep("sample", clinicalfiles, invert = TRUE, value = TRUE)

    if (length(clindatfile) > 1) {
        ncols <- vapply(clindatfile, function(x)
            ncol(readr::read_tsv(x, n_max = 3L)), integer(1L))
        clindatfile <- clindatfile[which.max(ncols)]
    }

    coldata <- cbioportal2clinicaldf(clindatfile)
    mdat <- cbioportal2metadata(mdatafile, licensefile)

    if (length(fusionExtra))
        fudat <- readr::read_tsv(fusionExtra, comment = "#")
    else
        fudat <- list()

    mdat <- c(mdat, metadats, fudat)

    if (any(.TCGAcols(coldata))) {
        gmap <- TCGAutils::generateMap(exptlist, coldata, TCGAbarcode)
    } else if (.hasMappers(coldata)) {
        gmap <- TCGAutils::generateMap(exptlist, coldata,
            sampleCol = "SAMPLE_ID", patientCol = "PATIENT_ID")
    } else {
        stop("Experiment data could not be mapped to colData")
    }

    MultiAssayExperiment(experiments = exptlist,
        colData = coldata, sampleMap = gmap, metadata = mdat)
}

cbioportal2metadata <- function(file, license) {
    md <- readLines(file, warn = FALSE)
    mdl <- lapply(seq_along(md), function(i) {
      sub(".+: ", "", md[[i]])
    })
    names(mdl) <- sub(":.+", "", md)
    lic <- readLines(license, warn = FALSE)
    lic <- paste0(lic[lic != ""], collapse = "\n")
    c(mdl, LICENSE = lic)
}

cbioportal2se <- function(file, ...) {
  df <- readr::read_tsv(file, comment = "#")
  looks.like.cn <- sapply(seq_along(df), function(i) {
    all(na.omit(df[[i]]) %in% -2:2)
  })
  numeric.cols <- sapply(df, class) == "numeric" | looks.like.cn
  rowdat <- DataFrame(df[,!numeric.cols])
  ##  rownames(rowdat) <- make.names(rowdat[, 1], unique=TRUE)
  se <-
    SummarizedExperiment(assays = as(df[, numeric.cols], "matrix"),
                                               rowData = rowdat)
  if(!all(grep("TCGA", rowData(se)[, 1])))
      return(NULL)
  rownames(se) <- rowData(se)[, 1]
  metadatafile <- sub("data", "meta", file)
  if(file.exists(metadatafile))
    metadata(se) <- cbioportal2metadata(metadatafile)
  return(se)
}

cbioportal2grl <-
  function(file,
           split.field,
           names.field) {
    df <- readr::read_tsv(file, comment = "#")
    if ("Strand" %in% colnames(df) && any(c(1L, -1L) %in% df$Strand)){
        newstrand <- rep("*", nrow(df))
        newstrand[df$Strand %in% 1L] <- "+"
        newstrand[df$Strand %in% -1L] <- "-"
        df$Strand <- newstrand
    }
    forbidden <-
      c(
        "seqnames",
        "ranges",
        "strand",
        "seqlevels",
        "seqlengths",
        "isCircular",
        "start",
        "end",
        "width",
        "element"
      )
    df <- df[,!colnames(df) %in% forbidden]
    split.field <-
      split.field[split.field %in% colnames(df)]
    if(grepl("\\.seg$", file))
        split.field <- colnames(df)[1]
    if (length(split.field) == 0)
      stop("No valid sample identifiers found")
    if (length(split.field) > 1) {
      stop(paste(
        "Select the correct identifier from:",
        paste(split.field, collapse = ", ")
      ))
    }
    names.field <- names.field[names.field %in% colnames(df)]
    if (length(names.field) == 0) {
      names.field <- NULL
    } else{
      names.field <- names.field[1]
    }
    chrom.col <- na.omit(match(c("Chromosome", "chrom"), colnames(df)))
    df <- df[!is.na(df[[chrom.col]]), ]
    grl <- makeGRangesListFromDataFrame(
      df,
      split.field = split.field,
      names.field = names.field,
      start.field = c("Start_Position", "loc.start", "Start"),
      end.field = c("End_Position", "loc.end", "End"),
      seqnames.field = c("Chromosome", "chrom"),
      strand.field = "Strand",
      keep.extra.columns = TRUE
    )
    if ("NCBI_Build" %in% colnames(df)) {
      genome(grl) <- df$NCBI_Build[1]
    }
    metadatafile <- sub("data", "meta", file)
    metadata(grl) <- cbioportal2metadata(metadatafile)
    grl <- RaggedExperiment::RaggedExperiment(grl)
    return(grl)
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
    rownames(clin) <- clin[["PATIENT_ID"]]
    return(clin)
}
