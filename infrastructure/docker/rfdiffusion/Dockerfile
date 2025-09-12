FROM nvcr.io/nvidia/cuda:11.6.2-cudnn8-runtime-ubuntu20.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get -q update \ 
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  git \
  wget \
  unzip \
  python3.9 \
  python3-pip \
  && python3.9 -m pip install -q -U --no-cache-dir pip \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get autoremove -y \
  && apt-get clean \
  && git clone https://github.com/RosettaCommons/RFdiffusion.git /app/RFdiffusion \
  && pushd /app/RFdiffusion \
  && git checkout 5606075d45bd23aa60785024b203ed6b0f6d2cd0 \
  && popd \
  && pip install -q --no-cache-dir \
  dgl==1.0.2+cu116 -f https://data.dgl.ai/wheels/cu116/repo.html \
  torch==1.12.1+cu116 --extra-index-url https://download.pytorch.org/whl/cu116 \
  e3nn==0.3.3 \
  pynvml==11.0.0 \
  decorator==5.1.0 \
  hydra-core==1.3.2 \
  pyrsistent==0.19.3 \
  icecream==2.1.3 \
  numpy==1.26.4 \
  /app/RFdiffusion/env/SE3Transformer \
  && wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
  && unzip awscliv2.zip \
  && ./aws/install \
  && rm awscliv2.zip

COPY run.sh setup.py /app/RFdiffusion/
WORKDIR /app/RFdiffusion
RUN pip install --no-cache-dir --no-deps . \
  && mkdir -p inputs outputs
ENV DGLBACKEND="pytorch"

ENTRYPOINT ["bash", "run.sh"]