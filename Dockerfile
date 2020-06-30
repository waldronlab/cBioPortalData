FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "cbioportaldata", "/bin/bash", "-c"]

COPY $PWD /app

RUN ./scripts/installdev.sh

