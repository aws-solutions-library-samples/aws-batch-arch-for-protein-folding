# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

ARG CUDA=11.3.1
ARG UBUNTU_VERSION=20.04
FROM nvcr.io/nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu${UBUNTU_VERSION}

# FROM directive resets ARGS, so we specify again
ARG CUDA

# Use bash to support string substitution.
SHELL ["/bin/bash", "-c"]

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
  cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
  libcusparse-dev-11-3 \
  libcublas-dev-11-3 \
  libcusolver-dev-11-3 \
  git \
  wget \
  unzip \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Install AWS CLI
RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m).zip" \
  && unzip awscliv2.zip \
  && ./aws/install \
  && rm awscliv2.zip

# Install Miniconda package manager and dependencies
ENV PATH="/opt/conda/bin:$PATH"
COPY environment.yml .

RUN wget -q -P /tmp -O /tmp/miniconda.sh \
  "https://repo.anaconda.com/miniconda/Miniconda3-py39_23.5.2-0-Linux-$(uname -m).sh" \
  && bash /tmp/miniconda.sh -b -p /opt/conda \
  && rm /tmp/miniconda.sh \
  && conda update -n base -c conda-forge conda \
  && conda config --set ssl_verify no \
  && conda init bash \
  && conda env update -n base --file environment.yml \
  && conda clean --all --yes

# Install OpenFold
RUN git clone https://github.com/aqlaboratory/openfold.git /opt/openfold \
  && cd /opt/openfold \
  # && git checkout 109442b14e6184fbee45e2696f21b052eb3fb1e5
  && git checkout 1d878a1203e6d662a209a95f71b90083d5fc079c

WORKDIR /opt/openfold
COPY setup.py .
RUN python setup.py install

# Install ESMFold
RUN git clone https://github.com/facebookresearch/esm.git /opt/esmfold \
  && cd /opt/esmfold \
  && git checkout 74d25cba46a7fd9a9f557ff41ed1d8e9f131aac3
WORKDIR /opt/esmfold
COPY run.sh /opt/esmfold/
COPY esmfold_inference.py /opt/esmfold/scripts/
RUN pip install .

RUN chmod 755 run.sh
ENTRYPOINT ["bash", "run.sh"]
