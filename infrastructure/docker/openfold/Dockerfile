# Original Copyright 2021 DeepMind Technologies Limited
# Modification Copyright 2022 # Copyright 2021 AlQuraishi Laboratory
# Modifications Copyright 2025 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM nvcr.io/nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04

# metainformation
LABEL org.opencontainers.image.version = "2.0.0"
LABEL org.opencontainers.image.authors = "OpenFold Team"
LABEL org.opencontainers.image.source = "https://github.com/aqlaboratory/openfold"
LABEL org.opencontainers.image.licenses = "Apache License 2.0"
LABEL org.opencontainers.image.base.name="docker.io/nvidia/cuda:12.4.1-devel-ubuntu22.04"

RUN apt-get update && apt-get install -y wget

RUN apt-key del 7fa2af80
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.0-1_all.deb
RUN dpkg -i cuda-keyring_1.0-1_all.deb

RUN apt-get install -y libxml2 cuda-minimal-build-12-1 libcusparse-dev-12-1 libcublas-dev-12-1 libcusolver-dev-12-1 git libaio-dev

RUN wget -P /tmp \
    "https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Miniforge3-Linux-x86_64.sh" \
    && bash /tmp/Miniforge3-Linux-x86_64.sh -b -p /opt/conda \
    && rm /tmp/Miniforge3-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

RUN git clone https://github.com/aqlaboratory/openfold.git /opt/openfold --depth=1
WORKDIR /opt/openfold

RUN mamba env update -n base --file /opt/openfold/environment.yml && mamba clean --all
RUN export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
RUN export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}
RUN wget -q -P /opt/openfold/openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

RUN python3 setup.py install

RUN git clone https://github.com/NVIDIA/cutlass --branch v3.6.0 --depth 1
RUN mamba env config vars set CUTLASS_PATH=$PWD/cutlass

# This setting is used to fix a worker assignment issue during data loading
RUN conda env config vars set KMP_AFFINITY=none

# COPY ./run_pretrained_openfold.py ./run.sh /opt/openfold/
COPY ./run.sh /opt/openfold/
ENTRYPOINT ["bash", "/opt/openfold/run.sh"]

