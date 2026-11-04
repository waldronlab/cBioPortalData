FROM ghcr.io/bioconductor/bioconductor_docker:devel

LABEL org.opencontainers.image.source="https://github.com/waldronlab/cBioPortalData"
LABEL org.opencontainers.image.description="Bioconductor devel with cBioPortalData installed"

# Dependencies only; this layer is cached until DESCRIPTION changes
COPY DESCRIPTION /tmp/pkg/DESCRIPTION
RUN Rscript -e "if (!requireNamespace('remotes', quietly = TRUE)) BiocManager::install('remotes', ask = FALSE, update = FALSE)" \
    && Rscript -e "deps <- remotes::local_package_deps('/tmp/pkg', dependencies = NA); BiocManager::install(c(deps, 'jsonlite', 'testthat'), ask = FALSE, update = FALSE)"

# The package itself
COPY . /tmp/pkg
RUN R CMD INSTALL --no-docs --no-multiarch /tmp/pkg \
    && rm -rf /tmp/pkg
