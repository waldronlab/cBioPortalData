# setwd("~/gh/cBioPortalData")
file_loc <- "inst/service/cBioPortal/api.json"

download.file(
    url = "https://www.cbioportal.org/api/v2/api-docs",
    destfile = file_loc
)

md5 <- digest::digest(file_loc, file = TRUE)
cbiolines <- readLines("R/cBioPortal.R")
mdline <- grep("\"[0-9a-f]{32}\"", cbiolines, value = TRUE)
oldmd5 <- unlist(strsplit(mdline, "\""))[[2]]
updatedlines <- gsub("\"[0-9a-f]{32}\"", dQuote(md5, FALSE), cbiolines)
github_output <- Sys.getenv("GITHUB_OUTPUT")

## success -- updated API files and MD5
if (!identical(oldmd5, md5)) {
    writeLines(updatedlines, con = file("R/cBioPortal.R"))
    if (nzchar(github_output)) {
        cat("api_changed=true\n", file = github_output, append = TRUE)
    }
    quit(status = 0)
} else {
    cat("::notice::API is unchanged; no update needed.\n")
    if (nzchar(github_output)) {
        cat("api_changed=false\n", file = github_output, append = TRUE)
    }
    quit(status = 0)
}
