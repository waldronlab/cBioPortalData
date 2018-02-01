#' @title Convert a data file downloaded from MSKCC's cBioPortal to
#' a MultiAssayExperiment object
#'
#'
#' @param cancer_study_id The cBioPortal study identifier
#' @param cancer_file The (optional) location of a previously downloaded tar
#' file
#' @param dir_location If no tar.gz file provided, the location to download
#' and untar the study files
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
#' mae <- importcBioPortal(laml, cancer_file = "data/laml_tcga.tar.gz")
#'
#' @export importcBioPortal
importcBioPortal <- function(cancer_study_id, cancer_file = NULL,
    dir_location = tempdir(), split.field = c("Tumor_Sample_Barcode", "ID"),
    names.field = c("Hugo_Symbol", "Entrez_Gene_Id")) {

    ## Load dataset to envir
    loc_data <- new.env(parent = emptyenv())
    data("studiesTable", envir = loc_data)
    studiesTable <- loc_data[["studiesTable"]]

    url_location <- paste0("https://media.githubusercontent.com/media/",
        "cBioPortal/datahub/master/public")

    if (!missing(cancer_study_id)) {
    if (!S4Vectors::isSingleString(cancer_study_id))
        stop("Provide a single study identifier")
    if (!cancer_study_id %in% studiesTable[["cancer_study_id"]])
        stop("Provide a valid study identifier")
    }

    if (is.null(cancer_file)) {
        cancer_file <- file.path(dir_location,
            paste0(cancer_study_id, ".tar.gz"))
        download.file(file.path(url_location, basename(cancer_file)),
            destfile = cancer_file)
    } else if (!file.exists(cancer_file)) {
        stop("'cancer_file' must exist")
    } else {
        dir_location <- dirname(cancer_file)
    }

    fileList <- untar(cancer_file, list = TRUE)
    datafiles <- grep(fileList, pattern = "data.+\\.(txt|seg)$", value = TRUE)
    datafiles <- c(datafiles, grep("meta_study", fileList, value = TRUE),
        grep("/LICENSE", fileList, value = TRUE))

    untar(cancer_file, files = datafiles, exdir = dir_location)

    exptfiles <- file.path(dir_location,
        grep("clinical|study|LICENSE", datafiles, invert = TRUE, value = TRUE))
    clinicalfiles <- file.path(dir_location,
        grep("clinical", datafiles, value = TRUE))
    mdatafile <- file.path(dir_location,
        grep("meta_study", datafiles, value = TRUE))
    licensefile <- file.path(dir_location,
        grep("/LICENSE", datafiles, value = TRUE))

    expnames <- sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))
    expseq <- seq_along(exptfiles)
    names(expseq) <- expnames

    exptlist <- lapply(expseq, function(i, files, xpnames) {
        fname <- files[[i]]
        message(paste0("Working on: ", fname))
        dat <- readr::read_delim(fname, delim = "\t")
        dat <- as(as.data.frame(dat, check.names = FALSE), "DataFrame")
        cexp <- xpnames[[i]]
        if (grepl("meth", cexp)) {
            .getMethyl(dat)
        } else if (grepl("gist", cexp)) {
            .getGISTIC(dat)
        } else {
            .biocExtract(dat)
        }
    }, files = exptfiles, xpnames = expnames)

    names(exptlist) <-
        sub(".*data_", "", sub("\\.txt", "", basename(exptfiles)))

    clindatfile <- grep("sample", clinicalfiles, invert = TRUE, value = TRUE)

    if (length(clindatfile) > 1) {
        ncols <- vapply(clindatfile, function(x)
            ncol(readr::read_tsv(x, n_max = 3L)), integer(1L))
        clindatfile <- clindatfile[which.max(ncols)]
    }

    coldata <- cbioportal2clinicaldf(clindatfile)
    mdat <- cbioportal2metadata(mdatafile, licensefile)
    gmap <- TCGAutils::generateMap(exptlist, coldata, TCGAbarcode)

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
  library(SummarizedExperiment)
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
    library(GenomicRanges)
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
