library(cBioPortalData)

cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

gps <- genePanels(cbioportal)[["genePanelId"]]

set.seed(100)
(tests <- sample(studies, 10))

cblist <- vector("list", length(tests))
debugonce(cBioPortalData)

for (stud in tests) {
    cblist[[stud]] <- tryCatch({
        cBioPortalData(cbioportal,
                studyId = stud, genePanelId = "IMPACT341")
    }, error = function(e) {
        conditionMessage(e)
    })
}


# save(studs, file = "liststudies.rda")
#
# (sum(sapply(studs, is, "try-error")) / length(studs)) * 100
#
# Filter(function(x) !is.null(x), lapply(studs, function(mae) {
#     if (is(mae, "MultiAssayExperiment")) {
#         list(length(mae), colnames(mae))
#     } else NULL
# }))
#
# "nccrcc_genentech_2014"
