# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM public.ecr.aws/lts/ubuntu:20.04

SHELL ["/bin/bash", "-c"]

# Install dependencies
RUN apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
      build-essential \
      cmake \
      wget \
      git \
      curl \
      tar \
      awscli \
      python3 \
      python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

# Install MMSeqs2 - Download pre-built static binary (x86_64 AVX2)
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz -O mmseqs.tar.gz \
  && tar xzf mmseqs.tar.gz \
  && mv mmseqs/bin/mmseqs /usr/local/bin/mmseqs \
  && chmod +x /usr/local/bin/mmseqs \
  && rm -rf mmseqs mmseqs.tar.gz

# Copy scripts
COPY . /opt/mmseqs

WORKDIR /opt/mmseqs

RUN chmod 755 /opt/mmseqs/run.sh \
  && chmod 755 /opt/mmseqs/run_mmseqs_msa.sh

ENTRYPOINT ["bash", "/opt/mmseqs/run.sh"]
