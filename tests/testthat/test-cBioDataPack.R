test_that("Check TGCA patient ID prefix removal works on non-char IDs", {
  df <- data.frame(PATIENT_ID=c(1,2,3,4))
  expect_equal(
    cBioPortalData:::.subBCLetters(df), df
  )
})

test_that("Check TGCA patient ID prefix removal works on character IDs", {
  df <- data.frame(PATIENT_ID=c('ABCDx00x0001', 'EFGHx00x0002'))
  expected <- data.frame(PATIENT_ID=c('TCGAx00x0001', 'TCGAx00x0002'))
  expect_equal(
    cBioPortalData:::.subBCLetters(df), expected
  )
})
