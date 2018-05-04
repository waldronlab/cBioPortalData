library(cgdsr)
library(S4Vectors)

cgds <- CGDS("http://www.cbioportal.org/")

clean_case_list <- function(study_id, cgds) {
    caseList <-
        getCaseLists(cgds, study_id)[, "case_list_description", drop = FALSE]
    newDesc <- tidyr::separate(caseList, "case_list_description",
        c("description", "samples"), sep = " \\((?=([0-9]))")
    newDesc[, "samples"] <- trimws(gsub("samples\\)", "", newDesc[["samples"]],
        ignore.case = TRUE))
    S4Vectors::DataFrame(cancer_study_id = study_id, newDesc)
}

## use coad to get all possible fields
coad <- clean_case_list("coadread_tcga_pub", cgds)

## copy and paste headers from http://www.cbioportal.org/data_sets.jsp
headers <- strsplit("All  Sequenced  CNA  RNA-Seq  Tumor mRNA (microarray)  Tumor miRNA  Methylation (HM27)  RPPA  Complete", "  ")[[1]]
coadnums <- strsplit("276  224  257  244  224  85  236  196  195", "  ")[[1L]]
names(coadnums) <- headers
coad$headers <- names(coadnums)[match(coad$samples, coadnums)]

headermap <- coad[, c("description", "headers")]
headermap[complete.cases(headermap), ]

prad <- clean_case_list("prad_tcga", cgds)

