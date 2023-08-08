#!/bin/bash

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# 

############################################################
# Deploy the AWS Batch Architecture for Protein Folding and Design in your AWS account
## Options
# -b S3 bucket name to use for deployment staging 
# -n CloudFormation stack name
# -r Deployment region
# -v ID of VPC to use. If left empty, a new VPC will be created.
# -s ID of private subnet to use. If left empty, a new VPC will be created.
# -d ID of default security group to use. If left empty, a new VPC will be created.
# -f ID of FSX for Lustre  to use. If left empty, a new FSx for Lustre instance will be created.
# -m Mount name of FSX for Lustre to use. If left empty, a new FSx for Lustre instance will be created.
# -i Create a SageMaker notebook instance? (Y/N)
# -e Automatically populate FSx for Lustre file system with model parameters and sequence databases? (Y/N)
# -g Create a compute environment for G5 instance types? (Y/N) Note that availabilty is region-specific.
# -a Use multiple availability zones? (Y/N)
# -p Create a compute environment for p4d.24xlarge instance types? (Y/N) Note that availabilty is region-specific.
# -c Name of the S3 bucket with the code repo zip file
#
# Example CMD
# ./deploy.sh \
#   -b "my-deployment-bucket" \
#   -n "my-batchfold-stack" \
#   -r "us-east-1" \
#   -e "Y" \
#   -g "Y" \
#   -c "my-deployment-bucket"

set -e
unset -v BUCKET_NAME STACK_NAME REGION VPC SUBNET \
    DEFAULT_SECURITY_GROUP FILE_SYSTEM_ID FILE_SYSTEM_MOUNT_NAME LAUNCH_SAGEMAKER_NOTEBOOK DOWNLOAD_FSX_DATA \
    CREATE_G5_COMPUTE_ENVIRONMENT MULTI_AZ CREATE_P4D_COMPUTE_ENVIRONMENT CODE_REPO_S3_BUCKET_NAME

TIMESTAMP=$(date +%s)

while getopts 'b:n:r:v:s:d:f:m:i:e:g:a:p:c:' OPTION; do
    case "$OPTION" in
    b) BUCKET_NAME="$OPTARG";;
    n) STACK_NAME="$OPTARG";;
    r) REGION="$OPTARG";;       
    v) VPC="$OPTARG";;        
    s) SUBNET="$OPTARG";;        
    d) DEFAULT_SECURITY_GROUP="$OPTARG";;        
    f) FILE_SYSTEM_ID="$OPTARG";;        
    m) FILE_SYSTEM_MOUNT_NAME="$OPTARG";;        
    i) LAUNCH_SAGEMAKER_NOTEBOOK="$OPTARG";;
    e) DOWNLOAD_FSX_DATA="$OPTARG";;        
    g) CREATE_G5_COMPUTE_ENVIRONMENT="$OPTARG";;
    a) MULTI_AZ="$OPTARG";;
    p) CREATE_P4D_COMPUTE_ENVIRONMENT="$OPTARG";;
    c) CODE_REPO_S3_BUCKET_NAME="$OPTARG";;
    *) exit 1 ;;    
    esac
done 

[ -z "$STACK_NAME" ] && { STACK_NAME="batchfold"; }
[ -z "$REGION" ] && { INPUT_FILE="us-east-1"; }
[ -z "$VPC" ] && { VPC=""; }
[ -z "$SUBNET" ] && { SUBNET=""; }
[ -z "$DEFAULT_SECURITY_GROUP" ] && { DEFAULT_SECURITY_GROUP=""; }
[ -z "$FILE_SYSTEM_ID" ] && { FILE_SYSTEM_ID=""; }
[ -z "$FILE_SYSTEM_MOUNT_NAME" ] && { FILE_SYSTEM_MOUNT_NAME=""; }
[ -z "$LAUNCH_SAGEMAKER_NOTEBOOK" ] && { LAUNCH_SAGEMAKER_NOTEBOOK="Y"; }
[ -z "$DOWNLOAD_FSX_DATA" ] && { DOWNLOAD_FSX_DATA="Y"; }
[ -z "$CREATE_G5_COMPUTE_ENVIRONMENT" ] && { CREATE_G5_COMPUTE_ENVIRONMENT="N"; }
[ -z "$MULTI_AZ" ] && { MULTI_AZ="Y"; }
[ -z "$CREATE_P4D_COMPUTE_ENVIRONMENT" ] && { CREATE_P4D_COMPUTE_ENVIRONMENT="N"; }
[ -z "$CODE_REPO_S3_BUCKET_NAME" ] && { CODE_REPO_S3_BUCKET_NAME="aws-batch-architecture-for-alphafold-public-artifacts"; }

zip -r code.zip * -x .\*/\* -x random_commands.sh -x analysis\* -x build\*
aws s3 cp code.zip s3://$BUCKET_NAME/main/batch-protein-folding.zip
rm code.zip

aws cloudformation package --template-file infrastructure/cloudformation/batch-protein-folding-cfn-root.yaml --output-template infrastructure/cloudformation/batch-protein-folding-cfn-packaged.yaml --region $REGION --s3-bucket $BUCKET_NAME --s3-prefix cfn 
aws cloudformation deploy --template-file infrastructure/cloudformation/batch-protein-folding-cfn-packaged.yaml --capabilities CAPABILITY_IAM --stack-name $STACK_NAME --region $REGION --parameter-overrides S3Bucket=$BUCKET_NAME \
  LaunchSageMakerNotebook=$LAUNCH_SAGEMAKER_NOTEBOOK VPC=$VPC Subnet=$SUBNET DefaultSecurityGroup=$DEFAULT_SECURITY_GROUP FileSystemId=$FILE_SYSTEM_ID FileSystemMountName=$FILE_SYSTEM_MOUNT_NAME \
  DownloadFsxData=$DOWNLOAD_FSX_DATA CreateG5ComputeEnvironment=$CREATE_G5_COMPUTE_ENVIRONMENT MultiAZ=$MULTI_AZ CreateP4dComputeEnvironment=$CREATE_P4D_COMPUTE_ENVIRONMENT CodeRepoS3BucketName=$CODE_REPO_S3_BUCKET_NAME Timestamp=$TIMESTAMP
rm infrastructure/cloudformation/batch-protein-folding-cfn-packaged.yaml