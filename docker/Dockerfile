# Image for running the R/greta pipeline on RunPod CPU pods.
#
# Matches the local greta 0.6 setup (doc/cv_run_plan.md, section 1): R 4.6.1,
# greta at 282944f, conda-forge python 3.12.14 with TensorFlow 2.21.0 and
# TensorFlow Probability 0.25.0 (docker/requirements.txt is the pip freeze of
# the local greta06-env). The repository is not baked in: clone it at a fixed commit on
# the pod and rsync data/ and temporary/ to /workspace. See docker/README.md.

FROM rocker/geospatial:4.6.1@sha256:e67737e30c66e6d2237570cd647e2bfec59315e533713473ab2ba391992bfa2c

# CRAN snapshot for the packages added here (posit package manager, binaries
# for Ubuntu noble)
ARG CRAN_DATE=2026-10-01
ARG NCPUS=2
ENV CRAN=https://p3m.dev/cran/__linux__/noble/${CRAN_DATE}

# system packages: OpenBLAS as R's BLAS and LAPACK (selected through Debian
# alternatives, so no LD_PRELOAD), sshd, and tools for the pod
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
       libopenblas0-pthread \
       openssh-server tmux rsync htop less curl ca-certificates \
  && update-alternatives --set libblas.so.3-x86_64-linux-gnu \
       /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 \
  && update-alternatives --set liblapack.so.3-x86_64-linux-gnu \
       /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3 \
  && rm -rf /var/lib/apt/lists/*

# single-threaded BLAS and OpenMP by default; set these on the pod to override.
# OpenBLAS reads them when R starts, so they must be in the environment, not
# Renviron
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

# python environment, as the local greta06-env: a conda-forge python (made with
# micromamba), and the pip packages of its pip freeze. Not a venv of the system
# python: under reticulate, setuptools' distutils shim does not load there, and
# tensorflow_probability fails to import
ARG MICROMAMBA_VERSION=2.9.0-0
COPY docker/requirements.txt /opt/requirements.txt
RUN curl -fsSL -o /usr/local/bin/micromamba \
       https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64 \
  && chmod +x /usr/local/bin/micromamba \
  && MAMBA_ROOT_PREFIX=/opt/micromamba micromamba create -y -p /opt/greta-env \
       -c conda-forge --override-channels \
       python=3.12.14 pip setuptools=84.0.0 wheel=0.48.0 packaging=26.3 \
  && /opt/greta-env/bin/python -m pip install --no-cache-dir \
       -r /opt/requirements.txt \
  && rm -rf /opt/micromamba /usr/local/bin/micromamba /opt/greta-env/conda-meta
# conda-meta is removed so that reticulate takes this as a plain python
# installation, and does not look for a conda binary to activate it with; and
# reticulate is stopped from falling back to a python of its own
ENV RETICULATE_PYTHON=/opt/greta-env/bin/python \
    RETICULATE_USE_MANAGED_VENV=no

# R packages, from R/packages.R and the library() and :: calls in the scripts
# the pipeline runs (fits, folds, predictions, figures, two-stage, scoring).
# Packages the base image already has are kept at its versions.
RUN echo "options(repos = c(CRAN = '${CRAN}'), Ncpus = ${NCPUS})" \
       >> "${R_HOME}/etc/Rprofile.site" \
  && MAKEFLAGS=-j${NCPUS} install2.r --error --skipinstalled --ncpus ${NCPUS} \
       -r "${CRAN}" \
       tidyverse readxl lme4 terra tidyterra tidygeocoder \
       future future.apply future.callr DHARMa Hmisc patchwork extraDistr \
       ggtext geodata sf TMB RcppEigen fmesher Matrix coda scales fields mgcv \
       ggrepel scoringRules posterior loo matrixStats RColorBrewer R6 \
       processx callr remotes reticulate tensorflow \
  && rm -rf /tmp/downloaded_packages

# greta at 282944f (0.6.0.9000); its dependencies from the snapshot
ARG GRETA_REF=282944f
RUN MAKEFLAGS=-j${NCPUS} Rscript -e "remotes::install_github( \
       'greta-dev/greta@${GRETA_REF}', upgrade = 'never', \
       dependencies = NA)" \
  && Rscript -e "stopifnot(packageVersion('greta') == '0.6.0.9000')" \
  && rm -rf /tmp/downloaded_packages

# sshd: key-only root login, keys from PUBLIC_KEY at start (docker/entrypoint.sh)
RUN mkdir -p /run/sshd /root/.ssh \
  && chmod 700 /root/.ssh \
  && rm -f /etc/ssh/ssh_host_* \
  && printf '%s\n' \
       'PermitRootLogin prohibit-password' \
       'PasswordAuthentication no' \
       'KbdInteractiveAuthentication no' \
       'PubkeyAuthentication yes' \
       > /etc/ssh/sshd_config.d/10-key-only.conf

# runpodctl, for a pod to remove itself when its job is done
ARG RUNPODCTL_VERSION=2.14.0
ARG RUNPODCTL_SHA256=2e0fd370a52a0fc7e43a6434a209348a4f6836fcdf1ad2b093b609e938138be9
RUN curl -fsSL -o /usr/local/bin/runpodctl \
       https://github.com/runpod/runpodctl/releases/download/v${RUNPODCTL_VERSION}/runpodctl-linux-amd64 \
  && echo "${RUNPODCTL_SHA256}  /usr/local/bin/runpodctl" | sha256sum -c - \
  && chmod +x /usr/local/bin/runpodctl

COPY docker/entrypoint.sh docker/run_pod_job.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh /usr/local/bin/run_pod_job.sh

WORKDIR /workspace
EXPOSE 22
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["sleep", "infinity"]
