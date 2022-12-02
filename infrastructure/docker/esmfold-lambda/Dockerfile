FROM public.ecr.aws/lambda/python:3.9 as build

RUN yum update -y \
    && yum install git gcc gcc-c++ python3-devel gzip tar wget make openssl-devel -y \
    && yum clean all \
    && rm -rf /var/cache/yum

RUN wget https://github.com/Kitware/CMake/releases/download/v3.25.1/cmake-3.25.1.tar.gz \
    && tar xzvf cmake-3.25.1.tar.gz \
    && cd cmake-3.25.1 \
    && ./bootstrap --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cmake --version

RUN python -m pip install -U --no-cache-dir pip wheel \
    && python -m pip install -U --no-cache-dir \
      --extra-index-url https://download.pytorch.org/whl/cpu torch

# Install CPU-Only version of OpenFold
RUN git clone https://github.com/aqlaboratory/openfold.git /opt/openfold \
    && cd /opt/openfold \
    && git checkout 109442b14e6184fbee45e2696f21b052eb3fb1e5
WORKDIR /opt/openfold
COPY openfold-cpu/setup.py /opt/openfold/
COPY openfold-cpu/amber_minimize.py /opt/openfold/openfold/np/relax/
COPY openfold-cpu/config.py /opt/openfold/openfold/
COPY openfold-cpu/primitives.py /opt/openfold/openfold/model/
COPY openfold-cpu/structure_module.py /opt/openfold/openfold/model/
RUN python -m pip install --no-cache-dir .

# Install ESMFold
RUN python -m pip install --no-cache-dir "fair-esm[esmfold]"

FROM public.ecr.aws/lambda/python:3.9
COPY --from=build /var/lang/lib/python3.9/site-packages ${LAMBDA_TASK_ROOT}

RUN mkdir -p /tmp/hub/checkpoints \
    && cd /tmp/hub/checkpoints \
    && curl -O "https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_3B_v1.pt" \
    && curl -O "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt" \
    && curl -O "https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t36_3B_UR50D-contact-regression.pt"

COPY app.py ${LAMBDA_TASK_ROOT}
WORKDIR ${LAMBDA_TASK_ROOT}
ENV TORCH_HOME="/tmp"
CMD ["app.lambda_handler"] 