devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

data_links <- vapply(full_ids, downloadcBioPortal, character(1L))

# remain <- full_ids[which(full_ids == "ov_tcga"):length(full_ids)]
# vapply(remain, downloadcBioPortal, character(1L))

dsize <- vapply(names(data_links), .checkSize, logical(1L))

redl <- names(which(!dsize))

stopifnot(!length(redl))

