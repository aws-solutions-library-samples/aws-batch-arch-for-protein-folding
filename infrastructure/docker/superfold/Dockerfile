FROM ubuntu:latest

# Set up the environment
ENV PATH=$PATH:/usr/local/cuda/bin

RUN apt-get -q update \
  && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
  wget \
  unzip \
  && wget --no-check-certificate -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
  && unzip awscliv2.zip \
  && ./aws/install \
  && rm awscliv2.zip

# Install dependencies
RUN apt-get update && \
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
RUN conda install -p /usr -y \
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
RUN pip install tensorflow-gpu==2.12

# Set environment variables for CUDA and cuDNN
ENV CUDA_HOME=/usr/local/cuda
ENV LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda/extras/CUPTI/lib64:$LD_LIBRARY_PATH
ENV PATH=/usr/local/cuda/bin:$PATH

# Download and set up superfold
RUN git clone https://github.com/rdkibler/superfold.git /opt/superfold && \
    cd /opt/superfold/mmalign && \
    g++ -static -O3 -ffast-math -lm -o MMalign MMalign.cpp && \
    cd /

# Set environment variable USER
ENV USER=ec2-user

# Copy run.sh to /opt/superfold
COPY info_collection.py /opt/superfold/
COPY run_superfold.py /opt/superfold/
COPY run.sh /opt/superfold/
RUN chmod +x /opt/superfold/run.sh
COPY config.json /opt/superfold/info_collection/config.json

# Clean up
RUN apt-get clean && \
    conda clean -a -y && \
    pip cache purge

# Run script
ENTRYPOINT ["bash", "/opt/superfold/run.sh"]
