devtools::install(".", dependencies = TRUE)

library(cBioPortalData)

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
studiesTable <- denv[["studiesTable"]]

studies <- studiesTable$cancer_study_id
studies <- stats::setNames(studies, studies)

complete <- vector("list", length(studies))
names(complete) <- studies

for (stud in studies) {
    message("Working on: ", stud)
    complete[[stud]] <- tryCatch({
        is(
            cBioDataPack(cancer_study_id = stud),
            "MultiAssayExperiment"
        )
    }, error = function(e) conditionMessage(e))
}

studiesTable$building <- unname(complete)
usethis::use_data(studiesTable, overwrite = TRUE)

