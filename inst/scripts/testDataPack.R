library(cBioPortalData)

cbioportal <- cBioPortal()
studies <- stats::setNames(nm = getStudies(cbioportal)[["studyId"]])

set.seed(1234)
# original sample in documentation
ff <- sample(studies, 100)

cacheLoc <- "~/data/cBioPortal"
if (!dir.exists(cacheLoc))
    dir.create(cacheLoc)

for (studyId in ff) {
    currFile <- file.path(cacheLoc, paste0(studyId, ".rda"))
    if (!file.exists(currFile)) {
        start <- proc.time()
        res <- tryCatch({
            cBioDataPack(cancer_study_id = studyId)
            }, error = function(e) conditionMessage(e))
        end <- proc.time()
        total <- end - start
        reslist <- list(res, total)
        assign(studyId, reslist)
        save(list = studyId, file = currFile)
    }
}

maestats <- lapply(studies, function(x) {
    mm <- tryCatch({
        get(x)
    }, error = function(e) conditionMessage(e))
    if (is(mm, "try-error") || is.character(mm)) {
        load(file.path(cacheLoc, paste0(x, ".rda")))
        mm <- get(x)
    }
    list(
        MAEclass = is(mm[[1L]], "MultiAssayExperiment"),
        time = mm[[2L]],
        lengths = length(mm[[1L]])
    )
})

maes <- sapply(maestats, `[[`, "MAEclass")
100 * sum(maes)/length(maes)
# success rate for cBioDataPack out of 273 studies
# 79.4

timings <- sapply(maestats, `[[`, "time")
summary(timings["elapsed", ])
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   2.420    6.972   24.250   78.546   97.242 1010.656


res <- vector("list", length(studies))
names(res) <- studies

datafiles <- file.path("~/data/cBioPortal", paste0(studies, ".rda"))
dataenv <- new.env(parent = emptyenv())

for (dts in datafiles[1:2]) {
    objname <- gsub(".rda", "", basename(dts))
    load(dts, dataenv)
    object <- dataenv[[objname]]
    res[[dts]] <- c(class = class(object[[1]]),
        length = length(object[[1]]), timing = round(object[[2]]["elapsed"],1)
    )
}
