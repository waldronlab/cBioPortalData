test_that("Check API is consistent with online resource", {
    skip_if_offline()

    api_url <- system.file(package = "cBioPortalData", "service",
        "cBioPortal", "api.yaml", mustWork = TRUE)
    tmp <- tempfile()
    download.file("https://www.cbioportal.org/api/api-docs", destfile = tmp)
    expect_equal(
        digest::digest(api_url, file = TRUE),
        digest::digest(tmp, file = TRUE)
    )
})