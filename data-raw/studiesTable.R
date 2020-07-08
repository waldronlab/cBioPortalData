## Script to create resource table with file names and links
library(cgdsr)
library(magrittr)
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

URL <- gsub("<\\/{0,1}[Aaibr]{1,2}>", "", studiesTable[["description"]]) %>%
    str_extract_all("<.*?>") %>% IRanges::CharacterList() %>%
    S4Vectors::endoapply(., function(x) {
    gsub("<[aA]\\s[Hh][Rr][Ee][Ff]=|\"|>", "", x)
     })
studiesTable[["description"]] <- gsub("<.*?>", "", studiesTable[["description"]])
studiesTable <- as(studiesTable, "DataFrame")
studiesTable[["URL"]] <- URL

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

usethis::use_data(studiesTable, overwrite = TRUE)
