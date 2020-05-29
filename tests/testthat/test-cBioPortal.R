test_that("digesting different API works", {
    expect_true(
        !identical(
            digest(list(cBioPortal())),
            digest(list(cBioPortal("beta")))
        )
    )
})
