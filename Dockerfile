FROM ghcr.io/bimberlabinternal/discvr-base:latest

ARG GH_PAT='NOT_SET'

ADD . /RIRA

RUN if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
    && cd /RIRA \
    && Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
    && R CMD build . \
    && R CMD INSTALL --build *.tar.gz \
    && rm -Rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& python3 -m celltypist.command_line --update-models --quiet
