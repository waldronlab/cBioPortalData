library(cBioPortalData)

denv <- new.env(parent = emptyenv())
# setwd("../../..")
load("./data/studiesTable.rda", envir = denv)
studiesTable <- denv[["studiesTable"]]


## PACK BUILD
message("PACK BUILD")
studies <- studiesTable$cancer_study_id
studies <- stats::setNames(studies, studies)

comp_pack <- vector("logical", length(studies))
names(comp_pack) <- studies

if (identical(system("hostname"), "supermicro") &&
    identical(Sys.getenv("USER"), "mramos")) {

    library(BiocParallel)
    registered()
    params <- MulticoreParam(
        workers = 60, stop.on.error = FALSE, progressbar = TRUE
    )

    comp_pack <- bplapply(X = studies, FUN = function(x) {
        if (identical(x, "ccrcc_utokyo_2013")) {
            FALSE
        } else {
            dats <- tryCatch({
                cBioPortalData::cBioDataPack(cancer_study_id = x, ask = FALSE)
            }, error = function(e) conditionMessage(e))
            is(dats, "MultiAssayExperiment")
        }
    }, BPPARAM = params)
    comp_pack <- unlist(comp_pack)

} else {

    for (pack_stud in studies) {
        message("Working on: ", pack_stud)
        ## avoid segfault
        else
            comp_pack[[pack_stud]] <- is(
                tryCatch({
                    cBioDataPack(cancer_study_id = pack_stud, ask = FALSE)
                }, error = function(e) conditionMessage(e)),
                "MultiAssayExperiment"
            )
    }

}

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
previous <- denv[["studiesTable"]]

if (!identical(previous[["pack_build"]], comp_pack)) {
    studiesTable[["pack_build"]] <- comp_pack
    usethis::use_data(studiesTable, overwrite = TRUE)
}

