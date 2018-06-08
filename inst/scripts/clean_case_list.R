library(S4Vectors)
library(tidyr)
library(cgdsr)

clean_case_list <- function(study_id, cgds) {
    caseList <-
        getCaseLists(cgds, study_id)[, "case_list_description", drop = FALSE]
    nsamps <- as.numeric(stringr::str_match(caseList[["case_list_description"]],
        "[0-9]{1,5}(?=\\ssamples)"))
    newDESC <- gsub(" \\(*[0-9]{1,5} samples\\)*|\\.$", "",
        caseList[["case_list_description"]])
    S4Vectors::DataFrame(cancer_study_id = study_id, description = newDESC,
        samples = nsamps)
}

