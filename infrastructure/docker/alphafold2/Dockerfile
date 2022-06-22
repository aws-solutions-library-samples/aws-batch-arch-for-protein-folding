# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

ARG CUDA=11.1.1
ARG ALPHAFOLD2_VERSION=v2.2.2
ARG UBUNTU_VERSION=20.04
# ARG AWS_DEFAULT_REGION
# ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI
FROM nvcr.io/nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu${UBUNTU_VERSION} as base_image

# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA
ARG ALPHAFOLD2_VERSION
# ARG AWS_DEFAULT_REGION
# ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI
### ---------------------------------------------

# Use bash to support string substitution.
SHELL ["/bin/bash", "-c"]

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      build-essential \
      cmake \
      cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \      
      git \
      hmmer \
      kalign \
      tzdata \
      wget \
      # Address https://ubuntu.com/security/CVE-2022-24407
      libsasl2-modules-sql=2.1.27+dfsg-2ubuntu0.1 \ 
      # Address https://ubuntu.com/security/CVE-2022-24407, https://ubuntu.com/security/CVE-2021-3711,
      # and https://ubuntu.com/security/CVE-2022-0778
      openssl=1.1.1f-1ubuntu2.13 \
      # Address https://ubuntu.com/security/CVE-2021-33910
      systemd=245.4-4ubuntu3.17 \
      # Address additional MEDIUM vulnerabilities
      apt dpkg e2fsprogs gcc libc6 gzip libgcrypt20 libzstd1 zstd lz4 libnettle7 slapd=2.4.49+dfsg-2ubuntu1.9 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install mambaforge package manager.
RUN wget -q -P /tmp \
  https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
  && bash /tmp/Mambaforge-Linux-x86_64.sh -b -p /opt/mamba \
  && rm /tmp/Mambaforge-Linux-x86_64.sh

# Install conda packages.
ENV PATH="/opt/mamba/bin:$PATH"
RUN mamba update -qy mamba \
    && mamba install -y -c conda-forge -c bioconda \
      openmm=7.5.1 \
      cudatoolkit=${CUDA_VERSION} \
      pdbfixer \
      pip \
      python=3.9 \
      boto3 \
      protobuf=3.20.1 \
      six=1.15.0 \
      bioconda::hhsuite=3.3.0 \
    && mamba clean --all

RUN git clone --branch ${ALPHAFOLD2_VERSION} --depth 1 https://github.com/deepmind/alphafold.git /app/alphafold

# Install pip packages.
RUN pip3 install --upgrade pip \
    && pip3 install -r /app/alphafold/requirements.txt \
    && pip3 install --upgrade \
      jax==0.2.14 \
      jaxlib==0.1.69+cuda$(cut -f1,2 -d. <<< ${CUDA} | sed 's/\.//g') \
      -f https://storage.googleapis.com/jax-releases/jax_releases.html \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Apply OpenMM patch.
WORKDIR /opt/mamba/lib/python3.9/site-packages
RUN patch -p0 < /app/alphafold/docker/openmm.patch

# Add SETUID bit to the ldconfig binary so that non-root users can run it.
RUN chmod u+s /sbin/ldconfig.real

# We need to run `ldconfig` first to ensure GPUs are visible, due to some quirk
# with Debian. See https://github.com/NVIDIA/nvidia-docker/issues/1399 for
# details.
# ENTRYPOINT does not support easily running multiple commands, so instead we
# write a shell script to wrap them up.
WORKDIR /app/alphafold
RUN echo $'#!/bin/bash\n\
ldconfig\n\
python /app/alphafold/run_alphafold.py "$@"' > /app/run_alphafold.sh \
  && chmod +x /app/run_alphafold.sh
ENTRYPOINT ["/app/run_alphafold.sh"]

#############################
# Start AWS-specific updates
#############################
FROM base_image
ARG AWS_DEFAULT_REGION
ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI

# Install the AWS CLI
RUN mamba install -y -c conda-forge awscli

# Copy stereo_chemical file from s3
RUN aws s3 cp \
  s3://aws-batch-architecture-for-alphafold-public-artifacts/stereo_chemical_props/stereo_chemical_props.txt \
  /app/alphafold/alphafold/common/

# Add batch job wrapper
COPY run_s3_batch_job.sh /app/run_s3_batch_job.sh

# Clean up
RUN rm -rf \
  /app/alphafold/docker \
  /app/alphafold/imgs \
  /app/alphafold/notebooks \
  /app/alphafold/scripts \
  && chmod +x /app/run_s3_batch_job.sh

# Set new entrypoint
ENTRYPOINT ["/app/run_s3_batch_job.sh"]