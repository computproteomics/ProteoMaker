FROM rocker/shiny:4.4.2
LABEL maintainer="Veit Schwaemmle <veits@bmb.sdu.dk>"
LABEL description="Docker image of VSClust implementation on top of shiny-server. The number of to-be-installed R packages requires patience when building this image."

RUN apt-get update && apt-get install -y r-cran-rjava  && apt-get clean
RUN apt-get install -y libglpk-dev && apt-get clean

RUN R -e "install.packages(c('BiocManager','remotes','devtools'), repos='http://cran.us.r-project.org'); \
  update.packages(ask=F);"

RUN R -e 'BiocManager::install(c("preprocessCore","qvalue","limma","SummarizedExperiment","cleaver", "S4Vectors"))'

RUN R -e 'install.packages("devtools");options(timeout=10000);devtools::install_github("rr-2/PeptideRanger")'
RUN R -e 'install.packages("devtools");devtools::install_github("computproteomics/polystest")'

# libraries for the shiny app
RUN R -e 'install.packages(c("rmarkdown","shinyWidgets","shinyBS","bslib","shinyalert","shinyjs","DT"))'

# Install ProteoMaker
RUN R -e 'install.packages("devtools");devtools::install_github("computproteomics/ProteoMaker")'

RUN rm -rf /srv/shiny-server
RUN mkdir /srv/shiny-server
COPY inst/shiny/*  /srv/shiny-server/
COPY inst/shiny/www /srv/shiny-server/www
RUN mkdir /srv/shiny-server/SimulatedDataSets; chmod a+w /srv/shiny-server/SimulatedDataSets
