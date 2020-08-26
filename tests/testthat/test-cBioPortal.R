test_that("Check API is consistent with online resource", {
    skip_if_offline()

    expect_silent(cBioPortal())
})
