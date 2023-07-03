test_that("cBioDataPack works on at least 70% of studies", {

    cbioportal <- cBioPortal()
    studies <- stats::setNames(nm = getStudies(cbioportal)[["studyId"]])

    isMAE <- structure(vector("logical", length(studies)), .Names = studies)

    for (stud in studies) {
        message("Working on: ", stud)
        result <- try({
            study <- cBioDataPack(cancer_study_id = stud, check_build = FALSE)
        })
        isMAE[stud] <- is(result, "MultiAssayExperiment")
        removePackCache(cancer_study_id = stud, dry.run = FALSE)
    }

    successrate <- (100 * sum(isMAE)) / length(isMAE)

    expect_true(successrate > 80)
})

test_that(".get_build_result is working", {
    local_mocked_bindings(
        `.loadReportData` = function(...) {
            list(api_build =
                data.frame(
                    studyId = c("abc", "def"), api_build = c(TRUE, FALSE)
                )
            )
        }
    )
    expect_true(
        .get_build_result("abc", "api_build")
    )
    expect_false(
        .get_build_result("def", "api_build")
    )
    expect_identical(
        .get_build_result("xyz", "api_build"), NA
    )

    expect_true(
        .is_study_id_building("abc", "api_build")
    )
    expect_true(
        .is_study_id_building("def", "api_build", ask = FALSE)
    )
    expect_true(
        .is_study_id_building("xyz", "api_build", ask = FALSE)
    )
})
