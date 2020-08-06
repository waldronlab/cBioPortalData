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
            is(
                cBioPortalData(
                    cbioportal, studyId = stud, genePanelId = "IMPACT341"
                ),
                "MultiAssayExperiment"
            )
        }, error = function(e) conditionMessage(e))
    }

    successrate <- (100 * sum(complete)) / length(complete)

    expect_true(successrate > 90)
})
