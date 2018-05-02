library(GenomicDataCommons)
library(TCGAutils)

que <- files() %>%
    filter( ~ cases.project.project_id == "TCGA-COAD" &
        data_category == "Copy Number Variation" &
        data_type == "Copy Number Segment")
q <- manifest(que)

test <- q$id[1:5]

ff <- gdcdata(test)

flist <- lapply(ff, function(filename) {
    readr::read_tsv(filename, comment = "#") })
ffs <- dplyr::bind_rows(flist)

bframe <- filenameToBarcode(paste0(unique(ffs$Sample), ".grch38.seg.txt"))
bcodes <- bframe[["aliquots.submitter_id"]]

ffs$Tumor_Sample_Barcode <-
    bframe$aliquots.submitter_id[
        match(paste0(ffs$Sample, ".grch38.seg.txt"), bframe$file_name)]

ragged <- MultiAssayExperimentData:::.biocExtract(ffs,
    c("Hugo_Symbol", "Entrez_Gene_Id", "Gene"))


gdc_clinical(test, include_list_cols = FALSE)

