## Script to create resource table with file names and links
library(cgdsr)
library(magrittr)
library(stringr)
library(tibble)
library(dplyr)

cgds <- CGDS("http://www.cbioportal.org/public-portal/")
studiesTable <- getCancerStudies(cgds)
cancer_file <- paste0(studiesTable[["cancer_study_id"]], ".tar.gz")
studiesTable <- add_column(studiesTable, cancer_file = cancer_file, .after = 1L)
studiesTable <- rename(studiesTable, study_name = name)

URL <- gsub("<\\/{0,1}[Aaibr]{1,2}>", "", studiesTable[["description"]]) %>%
    str_extract_all("<.*?>") %>% IRanges::CharacterList() %>%
    S4Vectors::endoapply(., function(x) {
    gsub("<[aA]\\s[Hh][Rr][Ee][Ff]=|\"|>", "", x)
     })
studiesTable[["description"]] <- gsub("<.*?>", "", studiesTable[["description"]])
studiesTable <- as(studiesTable, "DataFrame")
studiesTable[["URL"]] <- URL

devtools::use_data(studiesTable, overwrite = TRUE)
