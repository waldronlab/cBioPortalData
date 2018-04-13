library(MultiAssayExperimentData)
library(BiocParallel)
library(RTCGAToolbox)

setCache("/data/16tb/cbio")

data(studiesTable)
initSub <- studiesTable$cancer_study_id[1:5]

dlinks <- lapply(initSub, downloadcBioPortal)


