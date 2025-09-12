# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# ARG CUDA=11.1.1
ARG CUDA=11.6.0
ARG ALPHAFOLD2_VERSION=v2.3.2
FROM nvcr.io/nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu18.04
# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA
ARG ALPHAFOLD2_VERSION

# Use bash to support string substitution.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \ 
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  build-essential \
  cmake \
  cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
  git \
  hmmer \
  kalign \
  tzdata \
  wget \
  awscli \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get autoremove -y \
  && apt-get clean

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
  && mkdir /tmp/hh-suite/build \
  && pushd /tmp/hh-suite/build \
  && cmake -DHAVE_AVX2=1 -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
  && make -j 4 && make install \
  && ln -s /opt/hhsuite/bin/* /usr/bin \
  && popd \
  && rm -rf /tmp/hh-suite

# Install Miniconda package manager.
RUN wget -q -P /tmp \    
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
  && bash /tmp/Miniforge3-Linux-x86_64.sh -b -p /opt/conda \
  && rm /tmp/Miniforge3-Linux-x86_64.sh  

# Install conda packages.
ENV PATH="/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:$LD_LIBRARY_PATH"

RUN conda install --quiet --yes --channel conda-forge pip python=3.11 openmm=8.0.0 pdbfixer --override-channels \
  && conda install --quiet --yes --channel nvidia cuda=${CUDA_VERSION} --override-channels \
  && conda clean --all --force-pkgs-dirs --yes

# RUN git clone --branch ${ALPHAFOLD2_VERSION} --depth 1 https://github.com/deepmind/alphafold.git /app/alphafold
RUN git clone --depth 1 https://github.com/deepmind/alphafold.git /app/alphafold

RUN wget -q -P /app/alphafold/alphafold/common/ \
  https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Install pip packages.
RUN pip3 install --upgrade pip --no-cache-dir \
  && pip3 install -r /app/alphafold/requirements.txt --no-cache-dir \
  && pip3 install --upgrade --no-cache-dir \
  jax==0.4.26 \
  jaxlib==0.4.26+cuda12.cudnn89 \
  -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Add SETUID bit to the ldconfig binary so that non-root users can run it.
RUN chmod u+s /sbin/ldconfig.real

# Currently needed to avoid undefined_symbol error.
RUN ln -sf /usr/lib/x86_64-linux-gnu/libffi.so.7 /opt/conda/lib/libffi.so.7

# We need to run `ldconfig` first to ensure GPUs are visible, due to some quirk
# with Debian. See https://github.com/NVIDIA/nvidia-docker/issues/1399 for
# details.
# ENTRYPOINT does not support easily running multiple commands, so instead we
# write a shell script to wrap them up.
WORKDIR /app/alphafold
COPY run.sh /app/alphafold/run.sh
RUN echo $'#!/bin/bash\n\
  ldconfig\n\
  python /app/alphafold/run_alphafold.py "$@"' > /app/run_alphafold.sh \
  && chmod +x /app/run_alphafold.sh /app/alphafold/run.sh

ENTRYPOINT ["bash", "/app/alphafold/run.sh"]
