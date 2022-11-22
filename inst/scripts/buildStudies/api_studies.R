library(cBioPortalData)

## API BUILD
message("API BUILD")

cbioportal <- cBioPortal()
studies <- stats::setNames(nm = getStudies(cbioportal)[["studyId"]])

comp_api <- vector("logical", length(studies))
names(comp_api) <- studies

err_api <- vector("character", length(studies))
names(err_api) <- studies

if (identical(tolower(Sys.getenv("IS_BIOC_BUILD_MACHINE")), "true")) {

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
        ## clean up data
        removeDataCache(
            cbioportal, studyId = api_stud, genePanelId = "IMPACT341",
            dry.run = FALSE
        )
    }

} else if (identical(Sys.getenv("IS_SUPERMICRO_MACHINE", "TRUE"))) {

    library(BiocParallel)
    params <- MulticoreParam(
        workers = 64, stop.on.error = FALSE,
        progressbar = TRUE, jobname = "cBioPortal"
    )
    res_api <- bplapply(
        studies,
        function(api_stud) {
            message("Working on: ", api_stud)
            dats <- tryCatch({
                cBioPortalData(
                    cbioportal, studyId = api_stud, genePanelId = "IMPACT341"
                )
            }, error = function(e) conditionMessage(e))
            comp <- is(dats, "MultiAssayExperiment")
            if (!comp)
                err <- dats
            else
                err <- ""
            removeDataCache(
                cbioportal, studyId = api_stud, genePanelId = "IMPACT341",
                dry.run = FALSE
            )
            list(comp_api = comp, err_api = err)
        }, BPPARAM = params
    )
    comp_api <- vapply(res_api, `[[`, logical(1L), "comp_api")
    err_api <- vapply(res_api, `[[`, character(1L), "err_api")
}

err_api <- Filter(nzchar, err_api)
err_api_info <- lapply(setNames(nm = unique(err_api)),
    function(x) names(err_api)[err_api == x])
# table(err_api)
save(err_api_info, file = "inst/extdata/api/err_api_info.rda")

api_build <- rev(stack(comp_api))
names(api_build) <- c("studyId", "api_build")

denv <- new.env(parent = emptyenv())
api_file <- system.file(
    "extdata", "api", "api_build.rda",
    package = "cBioPortalData", mustWork = TRUE
)
load(api_file, envir = denv)
prev <- denv[["api_build"]]

if (!identical(prev, api_build)) {
    save(api_build, file = "inst/extdata/api/api_build.rda")
}
