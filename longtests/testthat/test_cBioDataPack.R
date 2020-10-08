test_that("cBioDataPack works on at least 70% of studies", {

    data(studiesTable)
    studies <- studiesTable$cancer_study_id
    studies <- stats::setNames(studies, studies)

    complete <- vector("list", length(studies))
    names(complete) <- studies

    for (stud in studies) {
        message("Working on: ", stud)
        complete[[stud]] <- tryCatch({
            cBioDataPack(cancer_study_id = stud)
        }, error = function(e) conditionMessage(e))
        removePackCache(cancer_study_id = stud, dry.run = FALSE)
    }

    isMAE <- vapply(
        complete, function(x) is(x, "MultiAssayExperiment"), logical(1L)
    )

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 70)
})
