## Emacs, make this -*- mode: sh; -*-

## bioc2u-based image: installs cBioPortalData and its Bioconductor/CRAN
## dependencies as pre-built Ubuntu .deb binaries via apt + bspm.
## Bioc2u extends r2u to cover all Bioconductor packages, so installation
## is fast and dependency resolution is handled by the system package manager.
##
## DISCLAIMER: This image is for TESTING PURPOSES ONLY. It combines the devel
## branch of cBioPortalData with a Bioconductor *release* base image, a
## combination that Bioconductor does not support. For regular use, pull the
## image built on bioconductor/bioconductor_docker:devel instead:
##   docker pull ghcr.io/waldronlab/cbioportaldata:devel
## (built from inst/docker/bioc/Dockerfile)
##
## Usage:
##   docker build -f inst/docker/bioc2u/Dockerfile -t cbioportaldata:bioc2u .
##   docker run --rm -it --security-opt seccomp=unconfined cbioportaldata:bioc2u

## Noble (Ubuntu 24.04) with r2u + bioc2u apt repositories and bspm pre-configured.
## Pin to a full tag for reproducibility, e.g. noble-bioc-3.21-r-4.5.0
ARG BIOC2U_TAG=noble
FROM ghcr.io/bioconductor/bioc2u-user:${BIOC2U_TAG}

LABEL org.opencontainers.image.source="https://github.com/waldronlab/cBioPortalData" \
      org.opencontainers.image.description="FOR TESTING ONLY: cBioPortalData (devel) on a bioc2u Bioconductor release base. For regular use, see ghcr.io/waldronlab/cbioportaldata:devel" \
      org.opencontainers.image.licenses="AGPL-3.0"

## ── System libraries ────────────────────────────────────────────────────────
## libssl / libcurl are needed by httr; libxml2 by xml2 / openssl chain.
## Installing as .deb means no source compilation.
RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

## ── R package dependencies (via apt / bspm) ─────────────────────────────────
## bspm::enable() is already active in the base image, so install.packages()
## transparently delegates to apt for any package that has a .deb binary.
## All heavy Bioconductor deps (MultiAssayExperiment, SummarizedExperiment,
## RaggedExperiment, GenomicRanges, BiocFileCache, TCGAutils, …) are resolved
## by the Ubuntu package manager — no compilation required.
RUN Rscript -e "\
    options(warn = 2); \
    pkgs <- c( \
        'BiocManager', \
        'MultiAssayExperiment', \
        'SummarizedExperiment', \
        'RaggedExperiment', \
        'GenomicRanges', \
        'GenomeInfoDb', \
        'IRanges', \
        'S4Vectors', \
        'BiocFileCache', \
        'TCGAutils', \
        'RTCGAToolbox', \
        'AnVIL', \
        'httr', \
        'readr', \
        'dplyr', \
        'tidyr', \
        'tibble', \
        'digest' \
    ); \
    install.packages(pkgs); \
"

## ── cBioPortalData itself ────────────────────────────────────────────────────
## Install the checked-out source (the build context) so the image reflects
## the code in this repo/PR. Any remaining dependencies declared in
## DESCRIPTION are resolved by remotes (and delegated to apt via bspm).
COPY . /tmp/cBioPortalData
RUN Rscript -e "install.packages('remotes'); \
        remotes::install_local('/tmp/cBioPortalData', dependencies = TRUE, upgrade = 'never')" \
    && rm -rf /tmp/cBioPortalData

## ── Smoke test ───────────────────────────────────────────────────────────────
RUN Rscript -e "library(cBioPortalData); cbio <- cBioPortal(); stopifnot(is(cbio, 'cBioPortal'))"

CMD ["R"]
