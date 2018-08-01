library(GenomicDataCommons)
library(TCGAutils)

.textConvert <- function(charvec) {
    paste0(capture.output(dput(charvec)), collapse = " ")
}

# catsel <- available_values("files", "data_category")
catsel <- c("Transcriptome Profiling", "Copy Number Variation",
    "Clinical", "Simple Nucleotide Variation")

datacats <- .textConvert(catsel)

# typesel <- available_values("files", "data_type")
typesel <- c("Gene Expression Quantification",
    "miRNA Expression Quantification",
    "Copy Number Segment", "Clinical Supplement",
    "Masked Somatic Mutation")

datatyps <- .textConvert(typesel)

# expsel <- available_values("files", "experimental_strategy")
expsel <- c("RNA-Seq", "miRNA-Seq", "Genotyping Array",
    "_missing", "WXS")

dataexp <- .textConvert(expsel)

que <- files() %>%
    filter(
        as.formula(paste0(
        "~ cases.project.project_id == \"TCGA-COAD\" &",
        " data_category %in% ", datacats, " &",
        " data_type %in% ", datatyps, " &",
        " experimental_strategy %in% ", dataexp
        ))
    )

q <- manifest(que)

test <- q$id[1:5]

ff <- gdcdata(test)





## Copy Number Segment && Copy Number Variation
qq <- files() %>%
    filter(~ cases.project.project_id == "TCGA-COAD" &
    data_category == "Copy Number Variation" &
    data_type == "Copy Number Segment")

q <- manifest(qq)

test <- q$id[1:5]

ff <- gdcdata(test)

flist <- lapply(ff, function(filename) {
    readr::read_tsv(filename, comment = "#") })
ffs <- dplyr::bind_rows(flist)

bframe <- UUIDtoBarcode(names(ff), "file_id", end_point = "center")
bcodes <- bframe[["cases.samples.portions.analytes.aliquots.submitter_id"]]

ffs$file_name <- rep(names(ff), lapply(flist, nrow))

ffs$Tumor_Sample_Barcode <- bcodes[match(ffs$file_name, bframe$file_id)]

ragged <- MultiAssayExperimentData:::.biocExtract(ffs, NULL)

gdc_clinical(test, include_list_cols = FALSE)

