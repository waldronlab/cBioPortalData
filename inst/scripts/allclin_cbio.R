library(cBioPortalData)

devtools::load_all()

allstudies <- .invoke_bind(
    cbioportal, "getAllStudiesUsingGET"
)

studyIds <- allstudies[["studyId"]]

pancans <- grepl("pan_can", studyIds, fixed = TRUE)

allIds <- studyIds[!pancans]
allIds <- setNames(allIds, allIds)

res <- lapply(allIds, function(study) {
    ptids <- .invoke_bind(cbioportal, "getAllPatientsInStudyUsingGET",
        studyId = study)[["patientId"]]
    ca <- .invoke_bind(cbioportal, "getAllClinicalAttributesInStudyUsingGET",
        studyId = study)
    ethrace <- c("ETHNICITY", "RACE")
    if (all(ethrace %in% ca[["clinicalAttributeId"]]))
        tidyr::spread(
            .invoke_bind(
                    cbioportal, "fetchAllClinicalDataInStudyUsingPOST",
                        clinicalDataType = "PATIENT",
                        studyId = study,
                        attributeIds = c("ETHNICITY", "RACE"),
                        ids = ptids,
                        projection = "SUMMARY"
                    ),
            clinicalAttributeId, value)
    else
        tibble::tibble()
})


results <- Filter(length, res)

## percentage of studies with RACE and ETHNICITY attributeIds
round((length(results)/length(allIds))*100, 1)

results

allptids <- unname(unlist(lapply(results, `[`, "patientId")))
dups <- duplicated(allptids)
duplist <- split(dups, rep(names(results), vapply(results, nrow, integer(1L))))
orderduplist <- duplist[names(results)]

stopifnot(identical(names(results), names(orderduplist)))

urest <- Map(function(x, y) {
    x[!y, ]
}, x = results, y = orderduplist)

urest <- Filter(nrow, urest)

allrows <- dplyr::bind_rows(urest)

length(unique(allrows[["studyId"]]))


race2 <- gsub("\\[|\\]", "", tolower(allrows[["RACE"]]))

amin <- grepl("american indian", race2, fixed = TRUE)
race2[amin] <- "american indian or alaska native"

cauc <- grepl("caucasian", race2, fixed = TRUE)
race2[cauc] <- "white"

asianos <- grepl("[^cauc]asian", race2)
race2 == "filipino"

race2[asianos] <- "asian"

hwpacific <- c("samoan", "hawaiian", "fiji islander")
race2[race2 %in% hwpacific] <- "native hawaiian or other pacific islander"

ab <- c("african american", "black")
race2[race2 %in% ab] <- "black or african american"

table(race2)
