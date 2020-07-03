library(cBioPortalData)
library(BiocFileCache)

# setCache("/data/16tb/cbio")
cbioportal <- cBioPortal()

(studies <- getStudies(cbioportal)[["studyId"]])

complete <- vector("list", length(studies))
names(complete) <- studies

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
load("inst/scripts/maelistfromcbio.rda")
# load("maelistfromcbio.rda")

chars <- vapply(complete, is.character, logical(1L))
complete[chars]

successrate <- (1 - length(complete[chars]) / length(complete) )* 100
successrate
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

