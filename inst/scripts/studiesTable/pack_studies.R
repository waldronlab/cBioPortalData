library(cBioPortalData)

denv <- new.env(parent = emptyenv())
load("./data/studiesTable.rda", envir = denv)
studiesTable <- denv[["studiesTable"]]


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


denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
previous <- denv[["studiesTable"]]

if (!identical(previous[["pack_build"]], comp_pack)) {
    studiesTable[["pack_build"]] <- comp_pack
    usethis::use_data(studiesTable, overwrite = TRUE)
}

