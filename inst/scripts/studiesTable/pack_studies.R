devtools::install(".", dependencies = TRUE)

library(cBioPortalData)

## PACK BUILD
message("PACK BUILD")
studies <- studiesTable$cancer_study_id
studies <- stats::setNames(studies, studies)

comp_pack <- vector("logical", length(studies))
names(comp_pack) <- studies

for (pack_stud in studies) {
    message("Working on: ", pack_stud)
    ## avoid segfault
    if (identical(pack_stud, "ccrcc_utokyo_2013"))
        comp_pack[[pack_stud]] <- FALSE
    else
        comp_pack[[pack_stud]] <- is(
            tryCatch({
                cBioDataPack(cancer_study_id = pack_stud, ask = FALSE)
            }, error = function(e) conditionMessage(e)),
            "MultiAssayExperiment"
        )
}

studiesTable[["pack_build"]] <- comp_pack

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
}

missingStudy <- studiesTable$cancer_study_id[
    !studiesTable$cancer_study_id %in% names(comp_api)
]

if (length(missingStudy))
    message("These datasets are not in the new API: ",
        paste0(missingStudy, collapse = ", "))

studiesTable[["api_build"]] <- comp_api[studiesTable$cancer_study_id]

