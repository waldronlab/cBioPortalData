# setwd("~/gh/cBioPortalData")
file_loc <- "inst/service/cBioPortal/api.json"
r_file <- "R/cBioPortal.R"
sha_regex <- "\"[0-9a-f]{64}\""

download.file(
    url = "https://www.cbioportal.org/api/v2/api-docs",
    destfile = file_loc
)

sha256 <- digest::digest(file_loc, algo = "sha256", file = TRUE)
cbiolines <- readLines(r_file)
shaline <- grep(sha_regex, cbiolines, value = TRUE)
if (length(shaline) != 1L)
    stop(
        "expected exactly one sha256 digest in ", r_file,
        ", found ", length(shaline)
    )
oldsha256 <- unlist(strsplit(shaline, "\""))[[2]]
updatedlines <- gsub(sha_regex, dQuote(sha256, FALSE), cbiolines)
github_output <- Sys.getenv("GITHUB_OUTPUT")

## success -- updated API files and sha256
if (!identical(oldsha256, sha256)) {
    writeLines(updatedlines, con = r_file)
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
