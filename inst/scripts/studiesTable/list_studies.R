devtools::install(".", dependencies = TRUE)

library(cBioPortalData)

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
studiesTable <- denv[["studiesTable"]]

studies <- studiesTable$cancer_study_id
studies <- stats::setNames(studies, studies)

complete <- vector("logical", length(studies))
names(complete) <- studies

for (stud in studies) {
    message("Working on: ", stud)
    ## avoid segfault
    if (identical(stud, "ccrcc_utokyo_2013"))
        complete[[stud]] <- FALSE
    else
        complete[[stud]] <- is(
            tryCatch({
                cBioDataPack(cancer_study_id = stud)
            }, error = function(e) conditionMessage(e)),
            "MultiAssayExperiment"
        )
}

studiesTable$building <- unname(complete)
usethis::use_data(studiesTable, overwrite = TRUE)

q("no", 0)
