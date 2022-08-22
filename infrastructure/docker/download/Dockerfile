# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM public.ecr.aws/amazonlinux/amazonlinux:latest
COPY scripts /
RUN amazon-linux-extras install epel -y \
  && yum update -y \
  && yum install aria2 tar rsync unzip -y \
  && curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
  && unzip awscliv2.zip \
  && ./aws/install
VOLUME /fsx

ENTRYPOINT ["usr/bin/bash"]