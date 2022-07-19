# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

FROM public.ecr.aws/amazonlinux/amazonlinux:latest
COPY scripts /
RUN amazon-linux-extras install epel -y \
  && yum update -y \
  && yum install aria2 tar rsync unzip awscli -y \
  && yum clean all \
  && rm -rf /var/cache/yum \
  && chmod 755 *.sh
VOLUME /fsx

ENTRYPOINT ["bash", "-c"]