devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

part_ids <- full_ids[full_ids != "blca_tcga"]

data_links <- vapply(part_ids, downloadcBioPortal, character(1L))



### Troubleshoot

downloadcBioPortal("blca_tcga")

