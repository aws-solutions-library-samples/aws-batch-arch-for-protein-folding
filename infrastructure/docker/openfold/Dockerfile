# Original Copyright 2021 DeepMind Technologies Limited
# Modification Copyright 2022 # Copyright 2021 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM nvcr.io/nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu18.04

RUN apt-key del 7fa2af80
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub

RUN apt-get update \ 
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \  
  wget \
  libxml2 \
  cuda-minimal-build-11-3 \
  libcusparse-dev-11-3 \
  libcublas-dev-11-3 \
  libcusolver-dev-11-3 \
  git \
  awscli \
  && rm -rf /var/lib/apt/lists/* \
  && apt-get autoremove -y \  
  && apt-get clean 

RUN wget -P /tmp \
    "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
    && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniconda3-latest-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH

RUN git clone https://github.com/aqlaboratory/openfold.git /opt/openfold \
  && cd /opt/openfold \
  && git checkout 109442b14e6184fbee45e2696f21b052eb3fb1e5

# installing into the base environment since the docker container wont do anything other than run openfold
RUN conda install -qy conda==4.13.0 \
  && conda env update -n base --file /opt/openfold/environment.yml \
  && conda clean --all --force-pkgs-dirs --yes

RUN wget -q -P /opt/openfold/openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN patch -p0 -d /opt/conda/lib/python3.7/site-packages/ < /opt/openfold/lib/openmm.patch

WORKDIR /opt/openfold
COPY . /opt/openfold/

RUN pip3 install --upgrade pip --no-cache-dir \
  && python3 setup.py install \
  && chmod 755 /opt/openfold/run.sh

ENTRYPOINT ["bash", "/opt/openfold/run.sh"]
