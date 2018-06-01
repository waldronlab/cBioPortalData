devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

data_links <- vapply(full_ids, downloadcBioPortal, character(1L))

source("studyTableMeta.R")

cs <- cgdsr::CGDS("http://www.cbioportal.org/")
names(full_ids) <- full_ids

allsamples <- lapply(full_ids[1:40], function(study_id) {
    headtbl <- headerMap(study_id, cs)
    headtbl[complete.cases(headtbl), c("cancer_study_id", "samples", "headers")]
})

allsamps <- do.call(rbind, allsamples)

allsamps <- as.data.frame(allsamps)

spread(allsamps, "headers", "samples")

