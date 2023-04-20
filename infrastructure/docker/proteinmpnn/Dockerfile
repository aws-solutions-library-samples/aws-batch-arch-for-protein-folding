FROM public.ecr.aws/amazonlinux/amazonlinux:2

RUN yum update -y \
    && yum install wget unzip git which -y \
    && yum clean all \
    && rm -rf /var/cache/yum

RUN amazon-linux-extras install python3.8 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1 \
    && python -m ensurepip --upgrade \
    && python -m pip install --upgrade --no-cache-dir \
      --extra-index-url https://download.pytorch.org/whl/cpu \
      torch numpy

RUN git clone https://github.com/dauparas/ProteinMPNN.git /opt/proteinmpnn \
    && cd /opt/proteinmpnn \
    && git checkout be1d37b6699dcd2283ab5b6fc8cc88774e2c80e9 \
    && rm -rf colab_notebooks examples inputs outputs training .git
    
RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m).zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm awscliv2.zip

WORKDIR /opt/proteinmpnn

COPY run.sh run.sh
COPY protein_mpnn_run.py /opt/proteinmpnn/

ENTRYPOINT ["bash", "run.sh"]
