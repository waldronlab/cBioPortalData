library(cBioPortalData)
cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

complete <- vector("logical", length(studies))
names(complete) <- studies

for (stud in studies) {
    message("Working on: ", stud)
    if (is.null(complete[[stud]])) {
        complete[[stud]] <- tryCatch({
            is(
                cBioPortalData(
                    cbioportal, studyId = stud, genePanelId = "IMPACT341"
                ),
                "MultiAssayExperiment"
            )
        }, error = function(e) conditionMessage(e))
    }
}

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
studiesTable <- denv[["studiesTable"]]

missingStudy <- studiesTable$cancer_study_id[!studiesTable$cancer_study_id %in% names(complete)]
if (length(missingStudy))
    message("These datasets are not in the new API: ", paste0(missingStudy, collapse = ", "))

studiesTable[, "api_build"] <- complete[studiesTable$cancer_study_id]
usethis::use_data(studiesTable, overwrite = TRUE)

q("no", 0)
