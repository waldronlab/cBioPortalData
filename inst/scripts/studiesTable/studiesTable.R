## Script to create resource table with file names and links
pkgs <- c("S4Vectors", "cgdsr", "stringr", "stringi", "tibble", "dplyr",
    "RCurl", "usethis")
ninst <- !vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)
if (any(ninst))
    BiocManager::install(pkgs[ninst])

library(S4Vectors)
library(cgdsr)
library(stringr)
library(stringi)
library(tibble)
library(dplyr)
library(RCurl)
library(usethis)

cgds <- CGDS("http://www.cbioportal.org/")
studiesTable <- getCancerStudies(cgds)
names(studiesTable) <-
    gsub("name", "study_name", names(studiesTable), fixed = TRUE)

parseURL <- gsub("<\\/{0,1}[Aaibr]{1,2}>", "", studiesTable[["description"]])
URL <- str_extract_all(parseURL, "<.*?>")
cleanURL <- lapply(URL, function(x) { gsub("<[aA]\\s[Hh][Rr][Ee][Ff]=|\"|>", "", x) })
cleanURL <- vapply(cleanURL, paste, character(1L), collapse = ", ")

studiesTable[["description"]] <- gsub("<.*?>", "", studiesTable[["description"]])
studiesTable <- as(studiesTable, "DataFrame")
studiesTable[["URL"]] <- cleanURL

fileURLs <- file.path("https://cbioportal-datahub.s3.amazonaws.com",
    paste0(studiesTable[["cancer_study_id"]], ".tar.gz"))

## Requires internet connection
validURLs <- vapply(fileURLs, RCurl::url.exists, logical(1L))

studiesTable <- studiesTable[validURLs, ]
changeCol <- vapply(studiesTable, function(x)
    any(stringi::stri_enc_mark(unlist(x)) == "native"), logical(1L))
if (any(changeCol))
    studiesTable[, changeCol] <-
        stringi::stri_enc_toascii(studiesTable[, changeCol])

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
oldStudies <- denv[["studiesTable"]]

errcode <- 1

if (!identical(studiesTable, oldStudies)) {
    errcode <- 0
    usethis::use_data(studiesTable, overwrite = TRUE)
}

q("no", errcode)
