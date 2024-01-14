FROM ghcr.io/bimberlabinternal/cellmembrane:latest

ARG GH_PAT='NOT_SET'

ENV RETICULATE_PYTHON=/usr/bin/python3

# NOTE: this is required when running as non-root. Setting MPLCONFIGDIR removes a similar warning.
ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp
ENV CELLTYPIST_FOLDER=/tmp

RUN pip3 install numba celltypist

# NOTE: this is also added to support running as non-root. celltypist needs to write in ~/
RUN mkdir /userHome && chmod -R 777 /userHome
ENV HOME=/userHome

ADD . /RIRA

RUN cd /RIRA \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    # Due to Matrix/SeuratObject: https://github.com/mojaveazure/seurat-object/issues/166
    && Rscript -e "install.packages('SeuratObject', ask = FALSE, force = TRUE, type = 'source', repos = 'https://cloud.r-project.org')" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
