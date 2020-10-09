library(cBioPortalData)

denv <- new.env(parent = emptyenv())
setwd("~/gh/cBioPortalData")
load("./data/studiesTable.rda", envir = denv)
studiesTable <- denv[["studiesTable"]]

## API BUILD
message("API BUILD")

cbioportal <- cBioPortal()
studies <- getStudies(cbioportal)[["studyId"]]

comp_api <- vector("logical", length(studies))
names(comp_api) <- studies

err_api <- vector("character", length(studies))
names(err_api) <- studies

for (api_stud in studies) {
    message("Working on: ", api_stud)
    dats <- tryCatch({
        cBioPortalData(
            cbioportal, studyId = api_stud, genePanelId = "IMPACT341"
        )
    }, error = function(e) conditionMessage(e))
    success <- is(dats, "MultiAssayExperiment")
    if (success)
        comp_api[[api_stud]] <- success
    else
        err_api[[api_stud]] <- dats
    ## try to free up memory
    gc()
}

save(err_api, file = "inst/extdata/err_api.rda")

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
prev <- previous[["api_build"]]

if (!identical(prev, api_comps)) {
    studiesTable[["api_build"]] <- api_comps
    usethis::use_data(studiesTable, overwrite = TRUE)
}
