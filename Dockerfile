FROM rocker/shiny:3.6.1

RUN apt-get update -y && apt-get install -y libssl-dev
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='https://cran.us.r-project.org')"
RUN R -e "BiocManager::install('multtest')"
RUN R -e "install.packages('Seurat', dependencies=TRUE, repos='https://cran.us.r-project.org')"

WORKDIR /srv/shiny-server/hsIsletVeresAetal

RUN chown -R shiny:shiny .

CMD ["/usr/bin/shiny-server.sh"]

