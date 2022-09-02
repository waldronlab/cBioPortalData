test_that("Check API is consistent with online resource", {
    skip_if_offline()

    expect_silent(cBioPortal())
})

test_that("cBioPortal integrity holds", {
    cbio <- cBioPortal()
    expect_true(
        is(cbio, "cBioPortal")
    )
    expect_gt(
        length(operations(cbio)), 50
    )
    expect_gt(
        length(searchOps(cbio, "get")), 30 
    )
})

test_that("getDataByGenes returns empty list when no data found", {
    cbio <- cBioPortal()
    muts <- getDataByGenes(
        api = cbio,
        studyId = "gbm_tcga_pub",
        genes = "ACTB",
        by = "hugoGeneSymbol",
        molecularProfileIds = "gbm_tcga_pub_mutations"
    )
    expect_identical(muts, structure(list(), names = character(0L)))

    mols <- getDataByGenes(
        api = cbio,
        studyId = "gbm_tcga_pub",
        genes = "ACTB",
        by = "hugoGeneSymbol",
        molecularProfileIds = "gbm_tcga_pub_cna_rae"
    )
    expect_identical(mols, structure(list(), names = character(0L)))
})

test_that("queryGeneTable returns a tibble structure", {
    cbio <- cBioPortal()
    feats <- queryGeneTable(api = cbio, by = "hugoGeneSymbol", genes = "ACTB")
    expect_identical(
        dim(feats), c(1L, 3L)
    )
    expect_identical(
        names(feats), c("entrezGeneId", "hugoGeneSymbol", "type")
    )
    expect_identical(
        feats[["entrezGeneId"]], 60L
    )
    expect_identical(
        feats[["hugoGeneSymbol"]], "ACTB"
    )
})

test_that("allSamples returns a tibble structure", {
    cbio <- cBioPortal()
    sampleIds <- allSamples(cbio, "gbm_tcga_pub")
    expect_true(
        tibble::is_tibble(sampleIds)
    )
    expect_true(
        "sampleId" %in% names(sampleIds)
    )
})