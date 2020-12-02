library(cBioPortalData)

denv <- new.env(parent = emptyenv())
# setwd("~/gh/cBioPortalData")
load("./data/studiesTable.rda", envir = denv)
studiesTable <- denv[["studiesTable"]]


## PACK BUILD
message("PACK BUILD")
studies <- studiesTable$cancer_study_id
studies <- stats::setNames(studies, studies)

comp_pack <- vector("logical", length(studies))
names(comp_pack) <- studies

err_pack <- vector("character", length(studies))
names(err_pack) <- studies

if (identical(system("hostname"), "supermicro") &&
    identical(Sys.getenv("USER"), "mramos")) {

    library(BiocParallel)
    registered()
    params <- MulticoreParam(
        workers = 30, stop.on.error = FALSE, progressbar = TRUE
    )

    res_pack <- bplapply(X = studies, FUN = function(x) {
        err <- character(1L)
        if (identical(x, "ccrcc_utokyo_2013")) {
            comp <- FALSE
        } else {
            dats <- tryCatch({
                cBioPortalData::cBioDataPack(cancer_study_id = x, ask = FALSE)
            }, error = function(e) conditionMessage(e))
            success <- is(dats, "MultiAssayExperiment")
            if (!success)
                err <- dats
            comp <- success
        }
        list(comp_pack = comp, err_pack = err)
    }, BPPARAM = params)
    comp_pack <- vapply(res_pack, `[[`, logical(1L), 1L)
    err_pack <- vapply(res_pack, `[[`, character(1L), 2L)

} else {

    for (pack_stud in studies) {
        message("Working on: ", pack_stud)
        ## avoid segfault
        if (identical(pack_stud, "ccrcc_utokyo_2013")) {
            comp_pack[[pack_stud]] <- FALSE
        } else {
            dats <- tryCatch({
                cBioDataPack(cancer_study_id = pack_stud, ask = FALSE)
            }, error = function(e) conditionMessage(e))
            success <- is(dats, "MultiAssayExperiment")
            if (success)
                comp_pack[[pack_stud]] <- success
            else
                err_pack[[pack_stud]] <- dats
        }
        ## try to free up memory
        gc()
        ## clean up data
        removePackCache(cancer_study_id = pack_stud, dry.run = FALSE)
    }
}

err_pack <- Filter(nchar, err_pack)
err_pack_info <- lapply(setNames(nm = unique(err_pack)),
    function(x) names(err_pack)[err_pack == x])
# table(err_pack)
save(err_pack_info, file = "inst/extdata/err_pack_info.rda")

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
previous <- denv[["studiesTable"]]
prev <- previous[["pack_build"]]

if (!identical(prev, comp_pack)) {
    studiesTable[["pack_build"]] <- comp_pack
    usethis::use_data(studiesTable, overwrite = TRUE)
}

