# setwd("~/gh/cBioPortalData")
file_loc <- "inst/service/cBioPortal/api.json"

download.file(
    url = "https://www.cbioportal.org/api/api-docs",
    destfile = file_loc
)

md5 <- digest::digest(file_loc, file = TRUE)
cbiolines <- readLines("R/cBioPortal.R")
mdline <- grep("\"[0-9a-f]{32}\"", cbiolines, value = TRUE)
oldmd5 <- unlist(strsplit(mdline, "\""))[[2]]
updatedlines <- gsub("\"[0-9a-f]{32}\"", dQuote(md5), cbiolines)

## success -- updated API files and MD5
if (!identical(oldmd5, md5)) {
    writeLines(updatedlines, con = file("R/cBioPortal.R"))
    quit(status = 0)
} else {
## failure -- API the same
    quit(status = 1)
}

