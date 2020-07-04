test_that("cBioPortal API is working with most studies", {
    cbio <- cBioPortal()
    studies <- getStudies(cbio)[["studyId"]]
    # may cause segfault: investigation pending
    studies <- studies[studies != "ccrcc_utokyo_2013"]

    complete <- vector("list", length(studies))
    names(complete) <- studies

    for (stud in studies) {
        message("Working on: ", stud)
        complete[[stud]] <- tryCatch({
            cBioPortalData(
                cbioportal, studyId = stud, genePanelId = "IMPACT341"
            )
        }, error = function(e) conditionMessage(e))
    }

    isMAE <- vapply(
        complete, function(x) is(x, "MultiAssayExperiment"), logical(1L)
    )

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 90)
})
