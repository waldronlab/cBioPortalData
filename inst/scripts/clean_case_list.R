library(S4Vectors)
library(tidyr)
library(cgdsr)

clean_case_list <- function(study_id, cgds) {
    caseList <-
        getCaseLists(cgds, study_id)[, "case_list_description", drop = FALSE]
    newDesc <- tidyr::separate(caseList, "case_list_description",
        c("description", "samples"), sep = " \\((?=([0-9]))")
    newDesc[, "samples"] <- trimws(gsub("samples\\)", "", newDesc[["samples"]],
        ignore.case = TRUE))
    S4Vectors::DataFrame(cancer_study_id = study_id, newDesc)
}

