test_that("cBioPortal API is working with most studies", {

    cbio <- cBioPortal()
    studies <- getStudies(cbio)[["studyId"]]

    isMAE <- structure(vector("logical", length(studies)), .Names = studies)

    for (api_stud in studies) {
        message("Working on: ", api_stud)
        result <- try({
            cBioPortalData(
                cbio, studyId = api_stud, genePanelId = "IMPACT341"
            )
        })
        isMAE[stud] <- is(result, "MultiAssayExperiment")
        removeDataCache(
            cbio, studyId = api_stud, genePanelId = "IMPACT341", dry.run = FALSE
        )
    }

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 80)
})
