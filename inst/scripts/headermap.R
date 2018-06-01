source("clean_case_list.R")
cgds <- cgdsr::CGDS("http://www.cbioportal.org/")
## use coad to get all possible fields
coad <- clean_case_list("coadread_tcga_pub", cgds)

## copy and paste headers from http://www.cbioportal.org/data_sets.jsp
headers <- strsplit("All  Sequenced  CNA  RNA-Seq  Tumor mRNA (microarray)  Tumor miRNA  Methylation (HM27)  RPPA  Complete", "  ")[[1]]
coadnums <- strsplit("276  224  257  244  224  85  236  196  195", "  ")[[1L]]
names(coadnums) <- headers
coad$headers <- names(coadnums)[match(coad$samples, coadnums)]

headermap <- coad[, c("description", "headers")]
headermap[nrow(headermap), "headers"] <- "CNA and Seq"
headermap <- headermap[complete.cases(headermap), ]
headermap <- rbind(headermap,
    DataFrame(description = c("All samples with methylation (HM450) data",
        "All tumor samples that have CNA and sequencing data"),
        headers = c("Methylation (HM450)", "CNA and Seq")))

headermap <- headermap[-c(which(headermap[["description"]] ==
    "All tumor samples that have mRNA, CNA and sequencing data"),
which(headermap[["description"]] ==
    "All samples with mRNA expression data" &
    headermap[["headers"]] == "Sequenced")), ]


# saveRDS(headermap, file = "../extdata/study_headers.rds")

