#!/bin/bash

export AWS_DEFAULT_REGION="us-west-2"
export ACCOUNT_ID="767398120215"
export IMAGE_REPO_NAME="superfold-local"
export IMAGE_TAG="latest"
export BUILD_CONTEXT="."

docker run -it --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v "$(pwd)":/LocalBuild/env \
  -e "IMAGE_NAME=aws/codebuild/standard:4.0" \
  amazon/aws-codebuild-local:latest \
  -i aws/codebuild/standard:4.0 \
  -a /LocalBuild/env \
  -s /LocalBuild/env
