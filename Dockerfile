FROM ghcr.io/bimberlabinternal/cellmembrane:latest

ADD . /RIRA

ENV RETICULATE_PYTHON=/usr/bin/python3

# NOTE: this is required when running as non-root. Setting MPLCONFIGDIR removes a similar warning.
ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp
ENV CELLTYPIST_FOLDER=/tmp

RUN pip3 install numba celltypist

# NOTE: this is also added to support running as non-root. celltypist needs to write in ~/
RUN mkdir /userHome && chmod -R 777 /userHome
ENV HOME=/userHome

RUN cd /RIRA \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
