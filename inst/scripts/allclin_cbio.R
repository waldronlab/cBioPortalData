library(cBioPortalData)

.invoke_bind <- cBioPortalData:::.invoke_bind

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

noData <- names(Filter(function(g) !length(g), res))
results <- Filter(length, res)

save(results, file ="allStudiesWithRE.rda")

## percentage of studies with RACE and ETHNICITY attributeIds
round((length(results)/length(allIds))*100, 1)

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
race2[amin] <- "Am. Indian / Hawaiian or P.I."

hwpacific <- c("samoan", "hawaiian", "fiji islander",
    "native hawaiian or other pacific islander")
race2[race2 %in% hwpacific] <- "Am. Indian / Hawaiian or P.I."

cauc <- grepl("caucasian", race2, fixed = TRUE)
cauc <- cauc | race2 == "white"
race2[cauc] <- "White"

asianos <- grepl("[^cauc]asian", race2)
asianos <- asianos | race2 %in% c("filipino", "laotian", "asian")

race2[asianos] <- "Asian"

ab <- c("african american", "black", "black or african american")
race2[race2 %in% ab] <- "Black or African Am."

unk <- c("not reported", "not evaluated")
race2[race2 %in% unk] <- "unknown"

other <- c("unknown", "other", "hispanic")
race2[race2 %in% other] <- "Other / Unknown"

table(race2)

allrows[["race"]] <- race2

library(ggplot2)

dat <- allrows %>% group_by(race) %>% summarize(counts = n())
dat <- mutate(dat, percentage = round(counts / sum(counts) * 100, 1))

dat_sort <- arrange(dat, percentage)
dat_sort$race <- factor(dat_sort$race, levels = dat_sort$race)
dat_sort <- dat_sort %>% mutate(database = "cBioPortal")

ggplot(data = dat_sort, aes(x = database, y = percentage, fill = race)) +
    geom_col() +
    geom_text(aes(label = paste0(percentage, "%")),
        position = position_stack(vjust = 0.5)) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal(base_size = 24) + theme(axis.text.x = element_blank()) +
    ylab("Percentage") +
    xlab(paste0(length(unique(allrows[["studyId"]])), " cBioPortal Studies")) +
    labs(fill = "Race")
