library(cBioPortalData)

denv <- new.env(parent = emptyenv())
load("./data/studiesTable.rda", envir = denv)
studiesTable <- denv[["studiesTable"]]

## API BUILD
message("API BUILD")
cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

comp_api <- vector("logical", length(studies))
names(comp_api) <- studies

for (api_stud in studies) {
    message("Working on: ", api_stud)
    comp_api[[api_stud]] <- is(
        tryCatch({
            cBioPortalData(
                cbioportal, studyId = api_stud, genePanelId = "IMPACT341"
            )
        }, error = function(e) conditionMessage(e)),
        "MultiAssayExperiment"
    )
    ## try to free up memory
    gc()
}

missingStudy <- studiesTable$cancer_study_id[
    !studiesTable$cancer_study_id %in% names(comp_api)
]

if (length(missingStudy))
    message("These datasets are not in the new API: ",
        paste0(missingStudy, collapse = ", "))

api_comps <- comp_api[studiesTable$cancer_study_id]

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
previous <- denv[["studiesTable"]]

if (!identical(previous[["api_build"]], api_comps)) {
    studiesTable[["api_build"]] <- api_comps
    usethis::use_data(studiesTable, overwrite = TRUE)
}
