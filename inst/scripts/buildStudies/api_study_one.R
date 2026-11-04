## Build a MultiAssayExperiment for a single cBioPortal study via the API
##
## Usage: Rscript api_study_one.R <studyId>
##
## Environment variables:
##   GENE_PANEL_ID        gene panel used to select genes (default IMPACT341)
##   RESULT_DIR           directory for the result JSON (default "results")
##   GITHUB_STEP_SUMMARY  set by GitHub Actions; markdown summary is appended
##
## Exits with status 1 when a MultiAssayExperiment is not created.

studyId <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(studyId) || !nzchar(studyId))
    stop("Provide a 'studyId' as the first argument")

genePanelId <- Sys.getenv("GENE_PANEL_ID", "IMPACT341")
resultDir <- Sys.getenv("RESULT_DIR", "results")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(cBioPortalData)
})

## profiles not indexed by gene are not part of cBioPortalData() output
notByGene <- c("STRUCTURAL_VARIANT", "GENERIC_ASSAY")

.mdEscape <- function(x) {
    x <- gsub("[\r\n]+", " ", x)
    gsub("|", "\\|", x, fixed = TRUE)
}

.dims <- function(mae, id) {
    tryCatch(
        paste(dim(mae[[id]]), collapse = " x "),
        error = function(e) NA_character_
    )
}

message("Working on: ", studyId)
cbio <- cBioPortal()

profs <- tryCatch(
    molecularProfiles(cbio, studyId = studyId),
    error = function(e) conditionMessage(e)
)
if (is.character(profs)) {
    message("Unable to list molecular profiles: ", profs)
    profs <- data.frame(
        molecularProfileId = character(0L),
        molecularAlterationType = character(0L),
        datatype = character(0L)
    )
}
profs <- as.data.frame(
    profs[, c("molecularProfileId", "molecularAlterationType", "datatype")]
)

start <- Sys.time()
mae <- tryCatch({
    cBioPortalData(
        api = cbio,
        studyId = studyId,
        genePanelId = genePanelId,
        check_build = FALSE,
        ask = FALSE
    )
}, error = function(e) conditionMessage(e))
seconds <- round(as.numeric(difftime(Sys.time(), start, units = "secs")))

ok <- is(mae, "MultiAssayExperiment")
error <- if (ok) "" else as.character(mae)
expNames <- if (ok) names(mae) else character(0L)

profs[["status"]] <- ifelse(
    profs[["molecularAlterationType"]] %in% notByGene, "skipped",
    ifelse(!ok, "failed",
        ifelse(profs[["molecularProfileId"]] %in% expNames,
            "included", "empty")
    )
)
profs[["dims"]] <- vapply(profs[["molecularProfileId"]], function(id) {
    if (id %in% expNames) .dims(mae, id) else NA_character_
}, character(1L), USE.NAMES = FALSE)

extra <- setdiff(expNames, profs[["molecularProfileId"]])

result <- list(
    studyId = studyId,
    genePanelId = genePanelId,
    api_build = ok,
    error = error,
    seconds = seconds,
    n_experiments = length(expNames),
    extra_experiments = extra,
    profiles = profs
)
jsonlite::write_json(
    result,
    path = file.path(resultDir, paste0(studyId, ".json")),
    auto_unbox = TRUE, pretty = TRUE, na = "null"
)

## markdown summary
icons <- c(
    included = "✅ included",
    empty = "⚪ empty for gene panel",
    skipped = "➖ skipped (not gene-indexed)",
    failed = "❌ build failed"
)
heading <- sprintf(
    "### %s `%s` (%s, %d s)",
    if (ok) "✅" else "❌", studyId, genePanelId, seconds
)
if (ok && !length(expNames))
    heading <- sprintf(
        "### ⚠️ `%s` (%s, %d s): no experiments",
        studyId, genePanelId, seconds
    )
md <- c(heading, "")
if (!ok)
    md <- c(md, paste0("**Error:** ", .mdEscape(error)), "")
if (nrow(profs)) {
    md <- c(md,
        "| molecularProfileId | alterationType | datatype | included | dims |",
        "|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s |",
            profs[["molecularProfileId"]],
            profs[["molecularAlterationType"]],
            profs[["datatype"]],
            icons[profs[["status"]]],
            ifelse(is.na(profs[["dims"]]), "–", profs[["dims"]])
        )
    )
} else {
    md <- c(md, "_No molecular profiles found._")
}
if (length(extra))
    md <- c(md, "",
        paste0(
            "Additional experiments: ",
            paste0(
                "`", extra, "` (",
                vapply(extra, .dims, character(1L), mae = mae), ")",
                collapse = ", "
            )
        )
    )
md <- c(md, "")

writeLines(md)
summaryFile <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summaryFile))
    cat(md, file = summaryFile, sep = "\n", append = TRUE)

if (!ok)
    quit(save = "no", status = 1L)
