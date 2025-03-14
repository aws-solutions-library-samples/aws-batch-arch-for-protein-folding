FROM nvcr.io/nvidia/cuda:12.6.3-runtime-amzn2023 as base

WORKDIR /home

RUN dnf update -y \
    && dnf install -y 'python3.11' \
    && dnf clean all \
    && rm -rf /var/cache/dnf/*

ENV VIRTUAL_ENV='/opt/venv'
ENV PATH="${VIRTUAL_ENV}/bin:$PATH"

##############################################################################
FROM base as build

ENV UV_COMPILE_BYTECODE=1
ENV UV_NO_CACHE=1
ENV UV_PYTHON_DOWNLOADS='never'
ENV UV_PYTHON=$PYTHON

COPY requirements.txt .

RUN dnf install -y tar gzip gcc 'python3.11-devel' \
    && curl -LsSf https://astral.sh/uv/install.sh | sh \
    && /root/.local/bin/uv venv $VIRTUAL_ENV \
    && /root/.local/bin/uv pip install --no-cache-dir -r requirements.txt

###############################################################################
FROM base as runtime

COPY --from=build $VIRTUAL_ENV $VIRTUAL_ENV
COPY scripts/ /home/scripts/
WORKDIR /home/scripts

ENV TMPDIR="/tmp"

COPY run.sh /home/scripts
RUN chmod 755 run.sh
ENTRYPOINT ["bash", "run.sh"]