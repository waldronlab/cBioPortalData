setwd("~/gh/cBioPortalData")
file_loc <- "inst/service/cBioPortal/api.json"

download.file(
    url = "https://www.cbioportal.org/api/api-docs",
    destfile = file_loc
)

md5 <- digest::digest(file_loc, file = TRUE)
cbiolines <- readLines("R/cBioPortal.R")
grep("\"[0-9a-f]{32}\"", cbiolines, value = TRUE)
md5
updatedlines <- gsub("\"[0-9a-f]{32}\"", dQuote(md5), cbiolines)

writeLines(updatedlines, con = file("R/cBioPortal.R"))
