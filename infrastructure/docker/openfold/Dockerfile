ARG CUDA=10.2
ARG UBUNTU_VERSION=18.04
# ARG AWS_DEFAULT_REGION
# ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI
FROM nvcr.io/nvidia/cuda:${CUDA}-cudnn8-runtime-ubuntu${UBUNTU_VERSION} as base_image

# FROM directive resets ARGS, so we specify again (the value is retained if
# previously set).
ARG CUDA
# ARG AWS_DEFAULT_REGION
# ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI

RUN apt-key del 7fa2af80
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub

# Use bash to support string substitution.
SHELL ["/bin/bash", "-c"]

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
      wget \
      git \
      cuda-command-line-tools-$(cut -f1,2 -d- <<< ${CUDA//./-}) \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install mambaforge package manager.
RUN wget -q -P /tmp \
  https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
  && bash /tmp/Mambaforge-Linux-x86_64.sh -b -p /opt/mamba \
  && rm /tmp/Mambaforge-Linux-x86_64.sh

RUN git clone --depth 1 https://github.com/aqlaboratory/openfold.git /opt/openfold

# Install mamba packages.
ENV PATH="/opt/mamba/bin:$PATH"
RUN mamba env update -n base --file /opt/openfold/environment.yml \
    && mamba clean --all

# Apply OpenMM patch.
RUN patch -p0 -d /opt/mamba/lib/python3.7/site-packages/ < /opt/openfold/lib/openmm.patch

WORKDIR /opt/openfold
RUN python3 setup.py install

#############################
# Start AWS-specific updates
#############################
FROM base_image
ARG AWS_DEFAULT_REGION
ARG AWS_CONTAINER_CREDENTIALS_RELATIVE_URI

# Install the AWS CLI
RUN mamba install -y -c conda-forge awscli \
    && mamba clean --all

# Copy stereo_chemical file from s3
RUN aws s3 cp \
  s3://aws-batch-architecture-for-alphafold-public-artifacts/stereo_chemical_props/stereo_chemical_props.txt \
  /opt/openfold/openfold/resources

# Add batch job wrapper
COPY run_s3_batch_job.sh /opt/run_s3_batch_job.sh

# Clean up
RUN rm -rf \
  /opt/openfold/imgs \
  /opt/openfold/tests \
  /opt/openfold/notebooks \
  && chmod +x /opt/run_s3_batch_job.sh

# Set new entrypoint
ENTRYPOINT ["/opt/run_s3_batch_job.sh"]