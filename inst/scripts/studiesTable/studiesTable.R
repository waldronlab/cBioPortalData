## the getStudies() functionality replaces the studiesTable

delayedAssign("studiesTable", value = {
    warning(
        "'studiesTable' dataset is defunct; see 'getStudies()'",
        call. = FALSE
    )
})

save("studiesTable", eval.promises=FALSE,
    file = "../../../data/studiesTable.rda")

