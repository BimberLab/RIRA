FROM ghcr.io/bimberlabinternal/cellmembrane:latest

ARG GH_PAT='NOT_SET'

ENV RETICULATE_PYTHON=/usr/bin/python3

# NOTE: this is required when running as non-root. Setting MPLCONFIGDIR removes a similar warning.
ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp
ENV CELLTYPIST_FOLDER=/tmp

# TODO: remove scikit-learn once this is resolved: https://github.com/Teichlab/celltypist/issues/75#issuecomment-1629757568
RUN pip3 install numba scikit-learn==1.2.2 celltypist

# NOTE: this is also added to support running as non-root. celltypist needs to write in ~/
RUN mkdir /userHome && chmod -R 777 /userHome
ENV HOME=/userHome

ADD . /RIRA

RUN cd /RIRA \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    # NOTE: Related to: https://github.com/satijalab/seurat/issues/7328. Should revert to a release once patched.
    # && Rscript -e "remotes::install_github('satijalab/seurat', ref='443ab86684253d9a7290c3d38c2bc1d8db021776');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
