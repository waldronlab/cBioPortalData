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

test_that(".biocExtract gives a RaggedExperiment for single-sample mutations", {
    maf <- S4Vectors::DataFrame(
        Hugo_Symbol = c("TP53", "KRAS"), NCBI_Build = "GRCh37",
        Chromosome = c("17", "12"), Start_Position = c(7577120L, 25398284L),
        End_Position = c(7577120L, 25398284L), Strand = "+",
        Variant_Classification = c("Missense_Mutation", "Nonsense_Mutation"),
        Tumor_Sample_Barcode = "S1"
    )
    re <- suppressMessages(
        cBioPortalData:::.biocExtract(maf, "Hugo_Symbol", character(0L))
    )
    expect_s4_class(re, "RaggedExperiment")
    expect_identical(rownames(re), c("TP53", "KRAS"))
    expect_identical(colnames(re), "S1")
    expect_identical(
        unname(RaggedExperiment::assay(re, "Variant_Classification")[, 1L]),
        c("Missense_Mutation", "Nonsense_Mutation")
    )
})
