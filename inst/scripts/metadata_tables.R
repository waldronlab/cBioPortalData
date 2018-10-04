devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

data_links <- vapply(full_ids, downloadStudy, character(1L))

source("studyTableMeta.R")

cs <- cgdsr::CGDS("http://www.cbioportal.org/")
names(full_ids) <- full_ids

allsamples <- lapply(full_ids, function(study_id) {
    headtbl <- headerMap(study_id, cs)
    headtbl[complete.cases(headtbl), c("cancer_study_id", "samples", "headers")]
})

allsamps <- do.call(rbind, allsamples)

allsamps <- as.data.frame(allsamps)

metadata_table <- spread(allsamps, "headers", "samples")

## check all ids have some metadata listing
all(full_ids %in% metadata_table[["cancer_study_id"]])

saveRDS(metadata_table, file = "../extdata/metadata_table.rds")

