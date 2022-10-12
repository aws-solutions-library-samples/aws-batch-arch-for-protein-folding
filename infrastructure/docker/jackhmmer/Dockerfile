# Original Copyright 2021 DeepMind Technologies Limited
# Modification Copyright 2022 # Copyright 2021 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM public.ecr.aws/lts/ubuntu:18.04_stable as base_image

SHELL ["/bin/bash", "-c"]

RUN apt-get update \
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
      build-essential \
      cmake \
      wget \
      git \
      unzip \
      hmmer \
      tar \
      awscli \
      python3.8 \
      python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \      
    && apt-get clean \
    && rm /usr/bin/python3 \
    && ln -s /usr/bin/python3.8 /usr/bin/python3

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
  && mkdir /tmp/hh-suite/build \
  && pushd /tmp/hh-suite/build \
  && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
  && make && make install \
  && ln -s /opt/hhsuite/bin/* /usr/bin \
  && popd
  # && rm -rf /tmp/hh-suite    

# Compile kalign2 from source
RUN wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz \
  && mkdir -p /tmp/kalign2/build \
  && tar -xvzf current.tar.gz -C /tmp/kalign2 \
  && pushd /tmp/kalign2 \
  && ./configure \
  && make && make install \
  && mv /usr/local/bin/kalign /usr/bin/kalign \
  && popd \
  && rm -rf /tmp/kalign2

COPY . /opt/msa
WORKDIR /opt/msa

RUN pip3 install --upgrade pip --no-cache-dir \
  && python3 -m pip install -r /opt/msa/requirements.txt --no-cache-dir \
  && chmod 755 /opt/msa/run.sh

ENTRYPOINT ["bash", "/opt/msa/run.sh"]