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
test_that("structuralVariantData returns data by molecular profile", {
    cbio <- cBioPortal()
    prof <- "brca_tcga_pan_can_atlas_2018_structural_variants"
    svs <- structuralVariantData(
        api = cbio, molecularProfileIds = prof, entrezGeneIds = 2064
    )
    expect_identical(names(svs), prof)
    expect_true(
        all(c("site1EntrezGeneId", "site2EntrezGeneId") %in% names(svs[[1L]]))
    )
    expect_true(
        all(
            svs[[1L]][["site1EntrezGeneId"]] == 2064L |
            svs[[1L]][["site2EntrezGeneId"]] == 2064L
        )
    )
    svsamp <- structuralVariantData(
        api = cbio, molecularProfileIds = prof, entrezGeneIds = 2064,
        sampleIds = c("TCGA-D8-A1X5-01", "TCGA-3C-AAAU-01")
    )
    expect_identical(unique(svsamp[[1L]][["sampleId"]]), "TCGA-D8-A1X5-01")
    expect_identical(
        structuralVariantData(
            api = cbio, molecularProfileIds = character(0L),
            entrezGeneIds = 2064
        ),
        structure(list(), names = character(0L))
    )
})

test_that("genericAssayData returns data by molecular profile", {
    cbio <- cBioPortal()
    prof <- "brca_tcga_pan_can_atlas_2018_armlevel_cna"
    samps <- c("TCGA-3C-AAAU-01", "TCGA-3C-AALI-01")
    meta <- genericAssayMeta(api = cbio, molecularProfileId = prof)
    expect_true("stableId" %in% names(meta))
    expect_true(all(c("1p_status", "1q_status") %in% meta[["stableId"]]))

    ga <- genericAssayData(
        api = cbio, molecularProfileIds = prof,
        genericAssayStableIds = c("1p_status", "1q_status"),
        sampleIds = samps
    )
    expect_identical(names(ga), prof)
    expect_true(
        all(ga[[1L]][["genericAssayStableId"]] %in% c("1p_status", "1q_status"))
    )
    expect_true(all(ga[[1L]][["sampleId"]] %in% samps))

    ## all stable IDs used when not specified
    gall <- genericAssayData(
        api = cbio, molecularProfileIds = prof, sampleIds = samps
    )
    expect_gt(
        length(unique(gall[[1L]][["genericAssayStableId"]])),
        length(unique(ga[[1L]][["genericAssayStableId"]]))
    )
})

test_that("fetchData routes profiles by molecularAlterationType", {
    cbio <- cBioPortal()
    st <- "brca_tcga_pan_can_atlas_2018"
    profs <- paste0(st, c("_structural_variants", "_armlevel_cna", "_gistic"))
    res <- fetchData(
        api = cbio, studyId = st, molecularProfileIds = profs,
        entrezGeneIds = 2064,
        sampleIds = c("TCGA-3C-AAAU-01", "TCGA-D8-A1X5-01"),
        genericAssayStableIds = "17q_status"
    )
    expect_true(all(c(profs[1L], profs[2L]) %in% names(res)))
    expect_true("site1EntrezGeneId" %in% names(res[[profs[1L]]]))
    expect_true("genericAssayStableId" %in% names(res[[profs[2L]]]))

    ## default excludes generic assay profiles
    alls <- suppressWarnings(fetchData(
        api = cbio, studyId = st, entrezGeneIds = 2064,
        sampleIds = c("TCGA-3C-AAAU-01", "TCGA-D8-A1X5-01")
    ))
    expect_false(any(names(alls) %in% c(profs[2L])))
    expect_true(profs[1L] %in% names(alls))
})
