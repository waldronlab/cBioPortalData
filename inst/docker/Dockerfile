FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY inst/docker/environment.yml .
RUN conda env update -n base --file environment.yml

# Make RUN commands use the new environment:
# SHELL ["conda", "run", "-n", "cbioportaldata", "/bin/bash", "-c"]

COPY . /app

RUN ["Rscript", "-e", "remotes::install_local(dependencies=FALSE, repos = BiocManager::repositories())"]
