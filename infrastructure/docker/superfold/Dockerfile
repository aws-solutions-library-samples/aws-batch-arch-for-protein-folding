FROM ubuntu:latest

# Set up the environment
ENV PATH=$PATH:/usr/local/cuda/bin

# Install dependencies
RUN apt-get update && \
    unzip \
    awscli \
    apt-get install -y libx11-6 libxau6 libxext6 g++ wget tar git build-essential && \
    rm /bin/sh && ln -s /bin/bash /bin/sh

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh && \
    bash /opt/miniconda.sh -b -u -p /usr && \
    conda update -n base -c defaults conda && \
    rm /opt/miniconda.sh && \
    rm -rf /usr/lib/terminfo && \
    conda install -p /usr -c conda-forge python=3.8 pip=20.2.4

# Install JAX with CUDA support
RUN pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Install additional dependencies via conda
RUN conda install -p /usr \
    -c pyg \
    -c pytorch \
    -c conda-forge \
    -c bioconda \
    -c nvidia \
    cudatoolkit \
    cudnn \
    cuda-nvcc && \
    conda install -p /usr \
    -c pyg \
    -c pytorch \
    -c conda-forge \
    -c bioconda \
    -c nvidia \
    attrs \
    biopython==1.81 \
    blas \
    blosc \
    bokeh \
    cryptography \
    dataclasses \
    decorator \
    distributed \
    h5py \
    idna \
    ipykernel \
    ipython \
    ipython_genutils \
    matplotlib \
    mock \
    more-itertools \
    numpy=1.23.5 \
    pandas \
    pandoc \
    pip \
    psutil \
    setuptools \
    tqdm \
    typing \
    typing-extensions \
    wheel \
    zict \
    python-dateutil \
    pyrsistent

# Install additional Python packages via pip
RUN pip install dm-haiku==0.0.4 absl-py dm-tree tensorflow ml-collections opt_einsum biopython==1.81

# Download and set up superfold
RUN git clone https://github.com/rdkibler/superfold.git /opt/superfold && \
    cd /opt/superfold/mmalign && \
    g++ -static -O3 -ffast-math -lm -o MMalign MMalign.cpp && \
    cd /

# Download and set up alphafold weights
#RUN mkdir -p /opt/weights /opt/weights/params && \
#    wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar -P /opt/weights/ && \
#    tar -xvf /opt/weights/alphafold_params_2022-12-06.tar -C /opt/weights/params/ && \
#    rm /opt/weights/alphafold_params_2022-12-06.tar && \
#    echo /opt/weights/ > /opt/superfold/alphafold_weights.pth

# Set environment variable USER
ENV USER=ec2-user

# Copy run.sh to /opt/superfold
COPY info_collection.py /opt/superfold
COPY run_superfold.py /opt/superfold
COPY run.sh /opt/superfold

# Clean up
RUN apt-get clean && \
    conda clean -a -y && \
    pip cache purge

# Run script
ENTRYPOINT ["bash", "opt/superfold/run.sh"]
