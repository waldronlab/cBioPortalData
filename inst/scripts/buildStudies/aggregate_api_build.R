## Aggregate per-study results from api_study_one.R
##
## Usage: Rscript aggregate_api_build.R
##
## Environment variables:
##   RESULT_DIR        directory of per-study result JSON files (default
##                     "results")
##   EXPECTED_STUDIES  JSON array of all study IDs dispatched; studies without
##                     a result file are counted as failures
##   PARTIAL_RUN       "true" when only a subset of studies was dispatched;
##                     the report is merged rather than replaced
##   UPDATE_REPORT     "true" to write inst/extdata/api/{api_build,
##                     err_api_info}.json
##   MIN_SUCCESS_RATE  percent of studies that must build (default 80)
##   GITHUB_*          set by GitHub Actions; used for job links, the step
##                     summary, and step outputs

resultDir <- Sys.getenv("RESULT_DIR", "results")
partial <- identical(tolower(Sys.getenv("PARTIAL_RUN")), "true")
update <- identical(tolower(Sys.getenv("UPDATE_REPORT")), "true")
minRate <- as.numeric(Sys.getenv("MIN_SUCCESS_RATE", "80"))

api_file_json <- "inst/extdata/api/api_build.json"
err_api_info_file <- "inst/extdata/api/err_api_info.json"

.mdEscape <- function(x) {
    x <- gsub("[\r\n]+", " ", x)
    x <- gsub("|", "\\|", x, fixed = TRUE)
    ifelse(nchar(x) > 200, paste0(substr(x, 1, 197), "..."), x)
}

## per-study results
files <- list.files(resultDir, pattern = "\\.json$", full.names = TRUE)
results <- lapply(files, jsonlite::fromJSON)
names(results) <- vapply(results, `[[`, character(1L), "studyId")

expected <- jsonlite::fromJSON(Sys.getenv("EXPECTED_STUDIES", "[]"))
if (!length(expected))
    expected <- names(results)
missing <- setdiff(expected, names(results))

.count <- function(res, status) {
    profs <- res[["profiles"]]
    if (!length(profs)) 0L else sum(profs[["status"]] %in% status)
}
summ <- data.frame(
    studyId = expected,
    api_build = vapply(expected, function(x) {
        isTRUE(results[[x]][["api_build"]])
    }, logical(1L)),
    error = vapply(expected, function(x) {
        if (x %in% missing) "No result (job cancelled or timed out)"
        else results[[x]][["error"]]
    }, character(1L)),
    seconds = vapply(expected, function(x) {
        if (x %in% missing) NA_real_ else as.numeric(results[[x]][["seconds"]])
    }, numeric(1L)),
    n_experiments = vapply(expected, function(x) {
        if (x %in% missing) 0L else as.integer(results[[x]][["n_experiments"]])
    }, integer(1L)),
    included = vapply(expected, function(x) {
        if (x %in% missing) 0L else .count(results[[x]], "included")
    }, integer(1L)),
    available = vapply(expected, function(x) {
        if (x %in% missing) 0L
        else .count(results[[x]], c("included", "empty", "failed"))
    }, integer(1L)),
    row.names = NULL
)

successrate <- round((100 * sum(summ[["api_build"]])) / nrow(summ), 1)
passed <- successrate > minRate

## compare with the committed report
prev <- jsonlite::fromJSON(api_file_json) |> as.data.frame()
prevStatus <- setNames(prev[["api_build"]], prev[["studyId"]])
known <- summ[["studyId"]] %in% names(prevStatus)
was <- prevStatus[summ[["studyId"]][known]]
now <- summ[["api_build"]][known]
broken <- summ[["studyId"]][known][was & !now]
fixed <- summ[["studyId"]][known][!was & now]
added <- summ[["studyId"]][!known]

## links to each study's job
.jobLinks <- function() {
    repo <- Sys.getenv("GITHUB_REPOSITORY")
    run <- Sys.getenv("GITHUB_RUN_ID")
    token <- Sys.getenv("GITHUB_TOKEN")
    if (!nzchar(repo) || !nzchar(run) || !nzchar(token))
        return(character(0L))
    api <- Sys.getenv("GITHUB_API_URL", "https://api.github.com")
    url <- sprintf(
        "%s/repos/%s/actions/runs/%s/attempts/%s/jobs",
        api, repo, run, Sys.getenv("GITHUB_RUN_ATTEMPT", "1")
    )
    links <- character(0L)
    page <- 1L
    repeat {
        resp <- httr::GET(url,
            query = list(per_page = 100L, page = page),
            httr::add_headers(
                Authorization = paste("Bearer", token),
                Accept = "application/vnd.github+json"
            )
        )
        httr::stop_for_status(resp)
        jobs <- httr::content(resp)[["jobs"]]
        if (!length(jobs))
            break
        nms <- vapply(jobs, `[[`, character(1L), "name")
        urls <- vapply(jobs, `[[`, character(1L), "html_url")
        links <- c(links, setNames(urls, sub(".* / ", "", nms)))
        page <- page + 1L
    }
    links
}
links <- tryCatch(.jobLinks(), error = function(e) {
    message("Unable to retrieve job links: ", conditionMessage(e))
    character(0L)
})

## markdown summary
.ids <- function(x) {
    if (length(x)) paste0("`", x, "`", collapse = ", ") else "none"
}
icon <- ifelse(!summ[["api_build"]], "❌",
    ifelse(summ[["n_experiments"]] == 0L, "⚠️", "✅"))
studyCol <- ifelse(
    summ[["studyId"]] %in% names(links),
    sprintf("[`%s`](%s)", summ[["studyId"]], links[summ[["studyId"]]]),
    sprintf("`%s`", summ[["studyId"]])
)
ord <- order(summ[["api_build"]], summ[["n_experiments"]] > 0L,
    summ[["studyId"]])
genePanelId <- if (length(results)) results[[1L]][["genePanelId"]] else ""

md <- c(
    sprintf("## %s cBioPortal API build: %s%% (%d / %d studies, %s)",
        if (passed) "✅" else "❌", successrate,
        sum(summ[["api_build"]]), nrow(summ), genePanelId),
    "",
    sprintf("Required success rate: > %s%%", minRate),
    "",
    paste0("**Newly broken:** ", .ids(broken)),
    "",
    paste0("**Newly fixed:** ", .ids(fixed)),
    "",
    paste0("**Not in previous report:** ", .ids(added)),
    "",
    "| | studyId | profiles included / available | seconds | error |",
    "|---|---|---|---|---|",
    sprintf("| %s | %s | %d / %d | %s | %s |",
        icon[ord], studyCol[ord],
        summ[["included"]][ord], summ[["available"]][ord],
        ifelse(is.na(summ[["seconds"]][ord]), "–",
            summ[["seconds"]][ord]),
        .mdEscape(summ[["error"]][ord])
    ),
    "",
    paste(
        "⚠️: a MultiAssayExperiment with no experiments.",
        "Available profiles exclude those not indexed by gene."
    ),
    ""
)

summaryFile <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summaryFile)) {
    cat(md, file = summaryFile, sep = "\n", append = TRUE)
} else {
    writeLines(md)
}

## update the build report
changed <- FALSE
if (update) {
    ran <- setdiff(expected, missing)
    report <- prev
    if (!partial)
        report <- report[report[["studyId"]] %in% expected, , drop = FALSE]
    report <- rbind(report,
        data.frame(
            studyId = setdiff(ran, report[["studyId"]]),
            api_build = logical(length(setdiff(ran, report[["studyId"]])))
        )
    )
    report[["api_build"]][match(ran, report[["studyId"]])] <-
        summ[["api_build"]][match(ran, summ[["studyId"]])]
    rownames(report) <- NULL

    errs <- jsonlite::fromJSON(err_api_info_file)
    errs <- setNames(
        rep(names(errs), lengths(errs)), unlist(errs, use.names = FALSE)
    )
    errs <- errs[names(errs) %in% report[["studyId"]] & !names(errs) %in% ran]
    failed <- summ[summ[["studyId"]] %in% ran & !summ[["api_build"]], ]
    errs <- c(errs, setNames(failed[["error"]], failed[["studyId"]]))
    err_api_info <- lapply(setNames(nm = unique(errs)),
        function(x) names(errs)[errs == x])

    if (!identical(prev, report)) {
        jsonlite::write_json(report, path = api_file_json, pretty = TRUE)
        changed <- TRUE
    }
    prevErrs <- jsonlite::fromJSON(err_api_info_file)
    if (!identical(prevErrs, err_api_info)) {
        jsonlite::write_json(err_api_info, path = err_api_info_file)
        changed <- TRUE
    }
}

outputFile <- Sys.getenv("GITHUB_OUTPUT")
if (nzchar(outputFile))
    cat(
        sprintf("success_rate=%s", successrate),
        sprintf("passed=%s", tolower(passed)),
        sprintf("changed=%s", tolower(changed)),
        file = outputFile, sep = "\n", append = TRUE
    )

message("Success rate: ", successrate, "%")
