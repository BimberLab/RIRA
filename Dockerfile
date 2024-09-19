FROM ghcr.io/bimberlab/discvr-base:latest

ARG GH_PAT='NOT_SET'

ADD . /RIRA

RUN cd /RIRA \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& R CMD build . \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds
