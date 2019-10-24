## A version has no caching
library(cBioPortalDataA)

cbioportal <- cBioPortal()

accA <- cBioPortalData(cbioportal, studyId = "acc_tcga",
    genePanelId = "IMPACT341")

save(accA, file = "accMethodA.rda")

gbmA <- cBioPortalData(cbioportal, studyId = "gbm_tcga",
    genePanelId = "IMPACT341")

save(gbmA, file = "gbmMethodA.rda")

unloadNamespace("cBioPortalDataA")

library(cBioPortalDataB)

cbioportal <- cBioPortal()

accB <- cBioPortalData(cbioportal, studyId = "acc_tcga",
genePanelId = "IMPACT341")

save(accB, file = "accMethodB.rda")

gbmB <- cBioPortalData(cbioportal, studyId = "gbm_tcga",
genePanelId = "IMPACT341")

save(gbmB, file = "gbmMethodB.rda")

unloadNamespace("cBioPortalDataB")

load("accMethodA.rda")
load("accMethodB.rda")
identical(accA, accB)
## TRUE

load("gbmMethodA.rda")
load("gbmMethodB.rda")
identical(gbmA, gbmB)
## TRUE


