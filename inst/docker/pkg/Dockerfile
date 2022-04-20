
FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest')); BiocManager::install('Bioconductor/BiocManager')"
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest')); BiocManager::install(update = TRUE, ask=FALSE)"
RUN Rscript -e "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/cran/__linux__/focal/latest')); devtools::install(dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

