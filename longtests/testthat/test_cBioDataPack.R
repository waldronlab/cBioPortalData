test_that("cBioDataPack works on at least 70% of studies", {

    cbioportal <- cBioPortal()
    studies <- stats::setNames(nm = getStudies(cbioportal)[["studyId"]])

    isMAE <- structure(vector("logical", length(studies)), .Names = studies)

    for (stud in studies) {
        message("Working on: ", stud)
        isMAE[[stud]] <- tryCatch({
            study <- cBioDataPack(cancer_study_id = stud)
            is(study, "MultiAssayExperiment")
        }, error = function(e) conditionMessage(e))
        removePackCache(cancer_study_id = stud, dry.run = FALSE)
    }

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 70)
})
