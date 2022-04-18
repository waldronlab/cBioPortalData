library(MultiAssayExperiment)
library(BiocFileCache)
devtools::load_all()

setCache("/data/16tb/cbio")

cbio <- cBioPortal()
studies <- stats::setNames(nm = getStudies(cbio)[["studyId"]])

data_links <- vapply(full_ids, downloadStudy, character(1L))

# remain <- full_ids[which(full_ids == "ov_tcga"):length(full_ids)]
# vapply(remain, downloadcBioPortal, character(1L))

dsize <- vapply(names(data_links), .checkSize, logical(1L))

(redl <- names(which(!dsize)))

if (length(redl))
    vapply(redl, downloadStudy, character(1L), force = TRUE)

alldata <- lapply(full_ids, function(x) tryCatch(importcBioPortal(x),
    error = function(e) e))

classes <- vapply(alldata, function(x) is(x, "MultiAssayExperiment"), logical(1L))
successrate <- (sum(classes)/length(classes)) * 100

isMAE <- vapply(full_ids, function(x) is(tryCatch(importcBioPortal(x),
    error = function(e) e), "MultiAssayExperiment"), logical(1L))
sr <- (sum(isMAE)/length(isMAE)) * 100

