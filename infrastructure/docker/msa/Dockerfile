FROM public.ecr.aws/lts/ubuntu:18.04_stable as base_image

SHELL ["/bin/bash", "-c"]

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      build-essential \
      cmake \
      wget \
      git \
      unzip \
      rsync \
      nano \
      hmmer \
      tar \
      aria2 \
      awscli \
      python3.7 \
    && apt-get clean \
    && ln -s /usr/bin/python3.7 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Compile HHsuite from source.
RUN git clone --branch v3.3.0 https://github.com/soedinglab/hh-suite.git /tmp/hh-suite \
    && mkdir /tmp/hh-suite/build \
    && pushd /tmp/hh-suite/build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/hhsuite .. \
    && make -j 4 && make install \
    && ln -s /opt/hhsuite/bin/* /usr/local/bin \
    && popd \
    && rm -rf /tmp/hh-suite    

# # Compile kalign2 from source
RUN wget http://msa.sbc.su.se/downloads/kalign/current.tar.gz \
  && mkdir -p /tmp/kalign2/build \
  && tar -xvzf current.tar.gz -C /tmp/kalign2 \
  && pushd /tmp/kalign2 \
  && ./configure \
  && make -j 4 && make install \
  && popd \
  && rm -rf /tmp/kalign2
  
# Add Tini (https://github.com/krallin/tini)
RUN wget -O "/tini" "https://github.com/krallin/tini/releases/download/v0.19.0/tini-arm64" \
  && chmod +x /tini

COPY . /opt/msa
WORKDIR /opt/msa

ENTRYPOINT ["/tini", "--", "bash", "-c"]