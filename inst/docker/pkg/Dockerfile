FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
ENV CRAN='https://p3m.dev/cran/__linux__/noble/latest'

RUN Rscript -e "BiocManager::install(update = TRUE, ask=FALSE)"
RUN Rscript -e "remotes::install_local(dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

