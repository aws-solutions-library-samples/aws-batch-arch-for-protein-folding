FROM public.ecr.aws/amazonlinux/amazonlinux:latest

RUN yum update -y \
    && yum groupinstall "Development Tools" -y \
    && yum install wget unzip -y \
    && yum clean all \
    && rm -rf /var/cache/yum

RUN amazon-linux-extras install python3.8

# # RUN pip3 install --upgrade pip wheel --no-cache-dir --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host files.pythonhosted.org

# RUN git clone --depth 1 https://github.com/dauparas/ProteinMPNN.git /opt/proteinmpnn
# RUN python -m ensurepip --upgrade
# RUN pip install --upgrade --no-cache-dir --extra-index-url https://download.pytorch.org/whl/cpu torch numpy
    
# RUN wget -O "awscliv2.zip" "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
#     && unzip awscliv2.zip \
#     && ./aws/install \
#     && rm awscliv2.zip

# # WORKDIR /opt/proteinmpnn

# # COPY run.sh run.sh

# # ENTRYPOINT ["bash", "run.sh"]
