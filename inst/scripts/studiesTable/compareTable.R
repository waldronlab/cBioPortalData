denv <- new.env(parent = emptyenv())
load("./data/studiesTable.rda", envir = denv)
latest <- denv[["studiesTable"]]

denv <- new.env(parent = emptyenv())
data("studiesTable", package = "cBioPortalData", envir = denv)
last <- denv[["studiesTable"]]

errcode <- 1

if (!identical(latest, last)) {
    errcode <- 0
    usethis::use_data(studiesTable, overwrite = TRUE)
}

q("no", errcode)
