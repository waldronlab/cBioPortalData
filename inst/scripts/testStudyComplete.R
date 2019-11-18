library(cBioPortalData)
library(BiocFileCache)

cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

gps <- genePanels(cbioportal)[["genePanelId"]]

set.seed(100)
(tests <- sample(studies, 200))

cblist <- vector("list", length(tests))
names(cblist) <- tests

for (stud in tests) {
    cblist[[stud]] <- tryCatch({
        cBioPortalData(cbioportal,
                studyId = stud, genePanelId = "IMPACT341")
    }, error = function(e) {
        conditionMessage(e)
    })
}

save(cblist, file = "maelistfromcbio.rda")
load("maelistfromcbio.rda")

chars <- vapply(cblist, is.character, logical(1L))
cblist[chars]

successrate <- (1 - length(cblist[chars]) / length(cblist) )* 100
## 62%


devtools::load_all()
debugonce(cBioPortalData)
debugonce(.portalExperiments)
debugonce(getDataByGenePanel)
debugonce(cBioPortalData:::.generateIdConvert)

