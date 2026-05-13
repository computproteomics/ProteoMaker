# syntax=docker/dockerfile:1.7
FROM rocker/shiny:4.5.2
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of VSClust implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

ENV CRAN_REPO=https://cloud.r-project.org

#RUN apt-get update && apt-get install -y r-cran-rjava  && apt-get clean
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    libglpk-dev \
    libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -q -e "options(repos = c(CRAN = Sys.getenv('CRAN_REPO')), timeout = 600, download.file.method = 'libcurl'); \
  install.packages(c('BiocManager', 'pak')); \
  update.packages(ask = FALSE)"

RUN R -q -e 'options(timeout = 600); BiocManager::install(c("preprocessCore", "qvalue", "limma", "SummarizedExperiment", "cleaver", "S4Vectors"), ask = FALSE, update = FALSE)'

#RUN R -e 'install.packages("pak");options(timeout=10000);pak::pak("rr-2/PeptideRanger")'
RUN --mount=type=secret,id=github_pat \
      if [ -f /run/secrets/github_pat ]; then \
        export GITHUB_PAT="$(cat /run/secrets/github_pat)"; \
      fi; \
      R -q -e "options(timeout = 600); pak::pak(c('rr-2/PeptideRanger', 'computproteomics/polystest'), upgrade = FALSE)"
#RUN R -e 'pak::pak("computproteomics/polystest")'

# libraries for the shiny app
RUN R -q -e 'options(repos = c(CRAN = Sys.getenv("CRAN_REPO")), timeout = 600, download.file.method = "libcurl"); install.packages(c("rmarkdown", "shinyWidgets", "shinyBS", "bslib", "shinyalert", "shinyjs", "DT"))'

# Install ProteoMaker
COPY . /tmp/ProteoMaker
RUN R -q -e "options(timeout = 600); pak::local_install('/tmp/ProteoMaker', upgrade = FALSE)"
#RUN R -q -e "pak::pak('computproteomics/ProteoMaker')"

RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
COPY inst/shiny/*  /srv/shiny-server/
COPY inst/shiny/www /srv/shiny-server/www
RUN mkdir /srv/shiny-server/SimulatedDataSets; chmod a+w /srv/shiny-server/SimulatedDataSets
