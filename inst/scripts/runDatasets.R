library(MultiAssayExperiment)
devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id
names(full_ids) <- full_ids

data_links <- vapply(full_ids, downloadcBioPortal, character(1L))

# remain <- full_ids[which(full_ids == "ov_tcga"):length(full_ids)]
# vapply(remain, downloadcBioPortal, character(1L))

dsize <- vapply(names(data_links), .checkSize, logical(1L))

(redl <- names(which(!dsize)))

if (length(redl))
    vapply(redl, downloadcBioPortal, character(1L), force = TRUE)

alldata <- lapply(full_ids[14:15], function(x) tryCatch(importcBioPortal(x),
    error = function(e) e))

