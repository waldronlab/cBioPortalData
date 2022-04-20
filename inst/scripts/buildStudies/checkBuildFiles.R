# setwd("~/gh/cBioPortalData")
apibuildfile <- file.info("inst/extdata/api/api_build.rda")
packbuildfile <- file.info("inst/extdata/pack/pack_build.rda")

library(lubridate)
api_recent <- Sys.Date() == date(ymd_hms(apibuildfile$mtime))
pack_recent <- Sys.Date() == date(ymd_hms(packbuildfile$mtime))

errcode <- if (any(c(api_recent, pack_recent))) 0 else 1
quit("no", errcode)
