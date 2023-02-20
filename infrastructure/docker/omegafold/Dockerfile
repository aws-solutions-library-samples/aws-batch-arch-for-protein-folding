ARG CUDA_VERSION=11.6.2
ARG UBUNTU_VERSION=20.04
FROM nvcr.io/nvidia/cuda:${CUDA_VERSION}-cudnn8-runtime-ubuntu${UBUNTU_VERSION}
ARG PYTHON_VERSION=3.10.7

# Use bash to support string substitution.
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        zlib1g-dev \
        libncurses5-dev \
        libgdbm-dev \
        libnss3-dev \
        libssl-dev \
        libreadline-dev \
        libffi-dev \
        wget \
        unzip \
        tar \
        git \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

RUN mkdir /python && cd /python \
    && wget https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz \
    && tar -xvf Python-${PYTHON_VERSION}.tgz \
    && cd Python-${PYTHON_VERSION} && ./configure --enable-optimizations \
    && make install \
    && cd / && rm -rf /python \
    && pip3 install --upgrade pip --no-cache-dir \
    && pip3 install --upgrade wheel --no-cache-dir

RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm awscliv2.zip

RUN git clone https://github.com/HeliXonProtein/OmegaFold.git /opt/OmegaFold \
    && cd /opt/OmegaFold \
    && git checkout 313c873ad190b64506a497c926649e15fcd88fcd

WORKDIR /opt/OmegaFold
COPY pipeline.py omegafold/pipeline.py
COPY run.sh run.sh
RUN chmod 755 run.sh \    
    && pip3 install --upgrade . --no-cache-dir

ENTRYPOINT ["bash", "run.sh"]