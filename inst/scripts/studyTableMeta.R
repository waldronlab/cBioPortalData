source("clean_case_list.R")

headerMap <- function(cancer_study_id, cgds_object) {

    headermap <- readRDS(system.file("extdata/study_headers.rds",
        package = "cBioPortalData", mustWork = TRUE))

    id_table <- clean_case_list(cancer_study_id, cgds_object)

    id_table$headers <- headermap$headers[
        match(id_table$description, headermap$description)]

    id_table <- id_table[!duplicated(id_table[, c("samples", "headers")]), ]
    id_table[!duplicated(id_table$headers), ]
}

