# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM public.ecr.aws/amazonlinux/amazonlinux:2
COPY scripts /scripts
RUN amazon-linux-extras install epel -y \
  && yum update -y \
  && yum install aria2 tar rsync unzip -y \
  && curl  "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m).zip" -o "awscliv2.zip" \
  && unzip awscliv2.zip \
  && ./aws/install \
  && chmod -R 755 /scripts

ENTRYPOINT ["bash", "-c"]