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
    && git checkout e61ecb7e3c32e630ff7a34d16c3a43fcf8f8a8bd \
    && rm -rf colab_notebooks examples inputs outputs training .git
    
RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m).zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm awscliv2.zip

WORKDIR /opt/proteinmpnn

COPY run.sh run.sh

ENTRYPOINT ["bash", "run.sh"]
