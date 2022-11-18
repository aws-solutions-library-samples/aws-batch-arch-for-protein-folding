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
    git \
    wget \
    unzip \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm awscliv2.zip

#clone the RF Design repo
RUN git clone --depth 1 https://github.com/RosettaCommons/RFDesign.git /opt/rfdesign \
  && mkdir /opt/rfdesign/input /opt/rfdesign/output

# Install Miniconda package manager.
RUN wget -q -P /tmp \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
  && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
  && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN conda update -n base -c conda-forge conda \
  && conda config --set ssl_verify no \
  && conda init bash \
  && conda clean -afy

# Note that we are using a different version of pytorch than recommended
# in RFDesign readme (suggested version is pytorch=1.10.1)
RUN conda install --name base \
    python=3.8 \
    pytorch=1.11 \
    dgl-cuda11.3 \
    pyg \
    numpy==1.21.5 scipy==1.7.3 requests==2.28.1 packaging pip \
    openmm==7.5.1 \
    pandas==1.5.1 \
    -c pytorch -c dglteam -c pyg -c conda-forge

RUN pip install --upgrade --no-cache-dir \
    https://github.com/openmm/pdbfixer/archive/refs/tags/v1.7.tar.gz \
    icecream==2.1.3 \
    lie_learn==0.0.1.post1 \
    opt_einsum==3.3.0 \
    e3nn==0.3.4 \
  && conda clean -afy

COPY run.sh /opt/rfdesign
COPY inpaint.py /opt/rfdesign/inpainting

WORKDIR /opt/rfdesign
RUN chmod 755 run.sh
ENTRYPOINT ["bash", "run.sh"]
