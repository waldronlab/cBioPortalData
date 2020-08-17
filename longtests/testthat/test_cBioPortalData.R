test_that("cBioPortal API is working with most studies", {
    cbio <- cBioPortal()
    studies <- getStudies(cbio)[["studyId"]]
    # may cause segfault: investigation pending
    studies <- studies[studies != "ccrcc_utokyo_2013"]

    complete <- vector("logical", length(studies))
    names(complete) <- studies

    for (api_stud in studies) {
        message("Working on: ", api_stud)
        complete[[api_stud]] <- is(
            tryCatch({
                cBioPortalData(
                    cbioportal, studyId = api_stud, genePanelId = "IMPACT341"
                )
            }, error = function(e) conditionMessage(e)),
            "MultiAssayExperiment"
        )
    }

    successrate <- (100 * sum(complete)) / length(complete)

    expect_true(successrate > 90)
})
