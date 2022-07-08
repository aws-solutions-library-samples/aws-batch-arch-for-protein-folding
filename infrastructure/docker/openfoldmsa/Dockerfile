FROM public.ecr.aws/lts/ubuntu:18.04_stable as base_image

SHELL ["/bin/bash", "-c"]

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      build-essential \
      cmake \
      wget \
      git \
      zip \
      nano \
      hmmer \
      tar \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && pushd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make -j 4 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/bin \
    && popd \
    && rm -rf /tmp/hh-suite    

# Compile kalign2 from source
RUN wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz \
  && mkdir -p /tmp/kalign2/build \
  && tar -xvzf current.tar.gz -C /tmp/kalign2 \
  && pushd /tmp/kalign2 \
  && ./configure \
  && make -j 4 && make install \
  && ln -s /opt/kalign2/bin/* /usr/bin \
  && popd \
  && rm -rf /tmp/kalign2

# Install mambaforge package manager.
RUN wget -q -P /tmp \
    https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-aarch64.sh \
    && bash /tmp/Mambaforge-Linux-aarch64.sh -b -p /opt/mamba \
    && rm /tmp/Mambaforge-Linux-aarch64.sh

ENV PATH="/opt/mamba/bin:$PATH"
RUN mamba update -qy mamba \
    && mamba install -y -c conda-forge -c bioconda \
      python=3.7 \
      openmm=7.5.1 \
      numpy==1.21.2 \
      scipy==1.7.1 \
      setuptools \
      pip \
      pdbfixer \
      dm-tree \
      requests \
      tqdm \
      PyYAML \
      typing-extensions \
    && mamba clean --all

COPY . /opt/openfold
WORKDIR /opt/openfold
RUN pip install .

#######################################################################################
# AWS-specific updates
#######################################################################################
FROM base_image
ARG AWS_DEFAULT_REGION
ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI

RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" \
  && unzip awscliv2.zip \
  && ./aws/install \
  && rm awscliv2.zip
  
# Add Tini (https://github.com/krallin/tini)
ENV TINI_VERSION v0.19.0
RUN wget -O "/tini" "https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini-arm64"  \
  && chmod +x /tini

# Clean up
RUN rm -rf \
  /opt/openfold/imgs \
  /opt/openfold/tests \
  /opt/openfold/notebooks

ENTRYPOINT ["/tini", "--"]