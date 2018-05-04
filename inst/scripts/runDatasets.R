devtools::load_all()

setCache("/data/16tb/cbio")

data(studiesTable)
full_ids <- studiesTable$cancer_study_id

data_links <- vapply(full_ids, downloadcBioPortal, character(1L))

# remain <- full_ids[which(full_ids == "ov_tcga"):length(full_ids)]
# vapply(remain, downloadcBioPortal, character(1L))

checkSize <- function(cancer_study_id) {

    bfc <- .get_cache()
    query_id <- glob2rx(cancer_study_id)
    study_file <- bfcquery(bfc, query_id, "rname")$rpath

    URL <- paste0("http://download.cbioportal.org/", cancer_study_id, ".tar.gz")

    header <- httr::HEAD(URL)$headers
    header_bytes <- as.numeric(header$`content-length`)

    local_bytes <- file.size(study_file)

    message("url: ", header_bytes, " vs. local: ", local_bytes)

    identical(header_bytes, local_bytes)
}

dsize <- vapply(names(data_links), checkSize, logical(1L))

redl <- names(which(!dsize))
