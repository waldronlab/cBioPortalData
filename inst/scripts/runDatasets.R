devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

data_links <- vapply(full_ids, downloadcBioPortal, character(1L))

