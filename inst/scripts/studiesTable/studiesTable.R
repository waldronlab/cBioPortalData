## the getStudies() functionality replaces the studiesTable

delayedAssign("studiesTable", value = {
    warning("'studiesTable' dataset is deprecated; see 'getStudies()'")
})

save("studiesTable", eval.promises=FALSE,
    file = "../../../data/studiesTable.rda")

