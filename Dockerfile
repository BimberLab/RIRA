FROM ghcr.io/bimberlabinternal/discvr-base:latest

ARG GH_PAT='NOT_SET'

ADD . /RIRA

RUN cd /RIRA \
    && if [ "${GH_PAT}" != 'NOT_SET' ];then echo 'Setting GITHUB_PAT to: '${GH_PAT}; export GITHUB_PAT="${GH_PAT}";fi \
	&& R CMD build . \
	# NOTE: remove zarr command when this is resolved: https://github.com/laminlabs/lamindb/pull/2334/files
	&& python3 -m pip install "zarr>=2.16.0,<3.0.0a0" \
	&& Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& python3 -m celltypist.command_line --update-models --quiet
