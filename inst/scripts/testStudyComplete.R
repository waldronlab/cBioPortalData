library(cBioPortalData)
library(BiocFileCache)

cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

complete <- vector("list", length(studies))
names(complete) <- studies

for (stud in studies) {
    complete[[stud]] <- tryCatch({
        cBioPortalData(cbioportal,
                studyId = stud, genePanelId = "IMPACT341")
    }, error = function(e) {
        conditionMessage(e)
    })
}

save(complete, file = "maelistfromcbio.rda")
load("maelistfromcbio.rda")

chars <- vapply(complete, is.character, logical(1L))
complete[chars]

successrate <- (1 - length(complete[chars]) / length(complete) )* 100
successrate
## 92.5

