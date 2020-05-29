test_that("digesting different API works", {
    expect_true(
        !identical(
            digest::digest(list(cBioPortal())),
            digest::digest(list(cBioPortal("beta.cbioportal.org")))
        )
    )
})
