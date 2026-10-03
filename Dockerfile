# One image with everything brava-annotate needs: VEP 105 + LOFTEE (GRCh38), bcftools/htslib,
# SpliceAI (TensorFlow 2.15) and the BRaVa scripts. Built from the same pixi.lock as the
# non-container install, so both give identical tool versions.
#
#   docker build -t brava-annotate .                       # CPU
#   docker build --build-arg PIXI_ENV=gpu -t brava-annotate:gpu .
FROM ubuntu:24.04 AS build
ARG PIXI_ENV=default
ARG PIXI_VERSION=v0.51.0
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl git && \
    curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=${PIXI_VERSION} PIXI_HOME=/usr/local bash
WORKDIR /opt/brava
COPY pixi.toml pixi.lock ./
RUN pixi install --locked -e ${PIXI_ENV} && \
    rm -rf /root/.cache/rattler && \
    # bioperl's meta-package drags in compilers, Java, docs and headers that nothing runs at runtime
    cd .pixi/envs/${PIXI_ENV} && \
    rm -rf docs include share/doc share/man share/info x86_64-conda-linux-gnu libexec/gcc lib/gcc lib/jvm && \
    find . -name '*.a' -delete && find . -name __pycache__ -prune -exec rm -rf {} +
# LOFTEE for GRCh38, pinned to the commit the original vep105_loftee image cloned
RUN git clone https://github.com/populationgenomics/loftee_38.git /opt/loftee && \
    git -C /opt/loftee checkout c8fdde00e515148450416128d43fcf01f1ee6bb8 && \
    rm -rf /opt/loftee/.git

FROM ubuntu:24.04
ARG PIXI_ENV=default
COPY --from=build /opt/brava/.pixi/envs/${PIXI_ENV} /opt/brava/.pixi/envs/${PIXI_ENV}
COPY --from=build /opt/loftee /opt/loftee
COPY bin /opt/brava/bin
COPY SAIGE_annotations/scripts /opt/brava/SAIGE_annotations/scripts
COPY resources /opt/brava/resources
COPY lib /opt/brava/lib
# Set the environment directly (no entrypoint) so `apptainer exec` and `docker run` behave the same
ENV CONDA_PREFIX=/opt/brava/.pixi/envs/${PIXI_ENV} \
    PATH=/opt/brava/bin:/opt/brava/.pixi/envs/${PIXI_ENV}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    PERL5LIB=/opt/loftee:/opt/brava/lib/perl \
    LOFTEE_PATH=/opt/loftee \
    BRAVA_HOME=/opt/brava \
    LC_ALL=C.UTF-8
CMD ["brava-annotate", "--help"]
