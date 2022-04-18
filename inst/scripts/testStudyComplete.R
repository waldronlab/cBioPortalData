library(cBioPortalData)
library(BiocFileCache)

# setCache("/data/16tb/cbio")
cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

complete <- structure(vector("list", length(studies)), .Names = studies)

for (stud in studies) {
    message("Working on: ", stud)
    if (is.null(complete[[stud]])) {
        complete[[stud]] <- tryCatch({
            cBioPortalData(
                cbioportal, studyId = stud, genePanelId = "IMPACT341"
            )
        }, error = function(e) conditionMessage(e))
    }
}

# save(complete, file = "inst/scripts/maelistfromcbio.rda")
# load("inst/scripts/maelistfromcbio.rda")
# load("maelistfromcbio.rda")

isMAE <- vapply(
    complete, function(x) is(x, "MultiAssayExperiment"), logical(1L)
)

successrate <- (sum(isMAE) * 100) / length(isMAE)
successrate
## success rate for cBioPortalData out of 268 studies
## 96.6

maelengths <- lengths(
    complete[
        vapply(complete, function(x) is(x, "MultiAssayExperiment"), logical(1L))
    ]
)

# save(maelengths, file = "maelenghtsfromcbio.rda")
# load("maelengthsfromcbio.rda")

table(maelengths)
# maelengths
#   1   2   3   4   5   6   7   8   9  10  11  12 13  14
# 103  47  10  11  38   5  10  20   4   2   4   2  2   1

