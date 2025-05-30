# Original Copyright 2021 DeepMind Technologies Limited
# Modification Copyright 2022 # Copyright 2021 AlQuraishi Laboratory
# Modifications Copyright 2025 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM nvidia/cuda:11.3.1-cudnn8-devel-ubuntu18.04

# metainformation
LABEL org.opencontainers.image.version = "1.0.0"
LABEL org.opencontainers.image.authors = "Gustaf Ahdritz"
LABEL org.opencontainers.image.source = "https://github.com/aqlaboratory/openfold"
LABEL org.opencontainers.image.licenses = "Apache License 2.0"
LABEL org.opencontainers.image.base.name="docker.io/nvidia/cuda:10.2-cudnn8-runtime-ubuntu18.04"

RUN apt-key del 7fa2af80
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub

RUN apt-get update && apt-get install -y wget libxml2 cuda-minimal-build-11-3 libcusparse-dev-11-3 libcublas-dev-11-3 libcusolver-dev-11-3 git
RUN wget -P /tmp \
    "https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Miniforge3-Linux-x86_64.sh" \
    && bash /tmp/Miniforge3-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniforge3-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

RUN git clone https://github.com/aqlaboratory/openfold.git /opt/openfold --depth=1
WORKDIR /opt/openfold

RUN mamba env update -n base --file /opt/openfold/environment.yml && mamba clean --all
RUN export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}
RUN wget -q -P /opt/openfold/openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN python3 setup.py install

COPY ./run_pretrained_openfold.py ./run.sh /opt/openfold/
ENTRYPOINT ["bash", "/opt/openfold/run.sh"]

