#' Title importcBioPortal: convert a .tar.gz file downloaded from
#' http://www.cbioportal.org/data_sets.jsp to a MultiAssayExperiment object
#'
#' @param tgzfile Path to the .tar.gz file as downloaded from cBioPortal web
#' page.
#' @param split.field A character vector of possible column names for the column
#' that is used to identify samples in a mutations or copy number file.
#' @param names.field A character vector of possible column names for the column
#' that is used to label ranges from a mutations or copy number file.
#'
#' @return A MultiAssayExperiment object
#' @author Levi Waldron
#' @export importcBioPortal
#'
#' @examples
#'
#' data(studiesTable)
#'
#' lamltar <- studiesTable[["cancer_study_id"]][3L]
#' mae <- importcBioPortal(lamltar)
importcBioPortal <- function(cancer_study_id, cancer_file = NULL,
                             split.field = c("Tumor_Sample_Barcode", "ID"),
                             names.field = c("Hugo_Symbol", "Entrez_Gene_Id"))
{
    url_location <-
      "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/"
    tempDir <- tempdir()
    untar(cancer_file, exdir = tempDir)
    fullpaths <- untar(tgzfile, list=TRUE)
    fullpaths <- grep(fullpaths,
            pattern = "data.+\\.(txt|seg)$", val=TRUE)
    setwd(unique(dirname(fullpaths)))
    datafiles <- dir(pattern = "data")
    datafiles = grep("clinical", datafiles, invert = TRUE, value = TRUE)
    exptlist <- lapply(datafiles, function(file) {
        message(paste0("Working on: ", file))
        fun <- get(.getFUN(file))
        fun(file, split.field=split.field, names.field=names.field)
    })
    names(exptlist) <-
        sub(".*data_", "", sub("\\.txt", "", basename(datafiles)))
    exptlist <- exptlist[!is.null(exptlist)]
    clindatfile <- grep("clinical", dir(pattern = "data"), val=TRUE)
    clindatfile <- grep("sample", clindatfile, invert=TRUE, val=TRUE)
    if(length(clindatfile) > 1){
        res <- sapply(clindatfile, function(x) ncol(readr::read_tsv(x)))
        clindatfile <- clindatfile[which.max(res)]
    }
    pdat <- cbioportal2clinicaldf(clindatfile)
    mdat <- cbioportal2metadata("meta_study.txt")
    setwd(orig.dir)
    library(MultiAssayExperiment)
    mae <-
        MultiAssayExperiment(experiments = exptlist,
                             colData = pdat,
                             metadata = mdat)
    return(mae)
}

cbioportal2metadata <- function(file) {
  file <- grep(file, dir(), val=TRUE)
  md <- readLines(file, warn = FALSE)
  mdl <- lapply(seq_along(md), function(i) {
    sub(".+: ", "", md[[i]])
  })
  names(mdl) <- sub(":.+", "", md)
  return(mdl)
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
  rownames(clin) <- clin$SAMPLE_ID
  return(clin)
}

.getFUN <- function(file) {
  cn <- colnames(read.delim(file, comment.char = "#", nrow = 1))
  is.gr <-
    any(c(
      "Start_Position",
      "loc.start",
      "Chromosome",
      "chrom",
      "Strand"
    ) %in% cn)
  ifelse(is.gr, "cbioportal2grl", "cbioportal2se")
}
