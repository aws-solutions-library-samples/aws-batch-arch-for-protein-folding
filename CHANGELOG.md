# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.13.2] - 2025-06-10

### 1.13.2 Fixed

- Updated build artifact location

---

## [1.13.1] - 2025-05-28

### 1.13.1 Fixed

- Fix issue with Download container build timing out

---

## [1.13.0] - 2025-03-14

### 1.13.0 Update

- Remove CodeCommit dependency

### 1.13.0 Fixed

- Resolve various cfn-lint issues
- Update OpenFold dockerfile

---

## [1.12.7] - 2024-03-14

### 1.12.7 Fixed

- Updated S3 location for reference data
- Added missing script to download ESM parameters

---

## [1.12.6] - 2023-12-08

### 1.12.6 Changed

- Remove jax from notebook requirements

---

## [1.12.5] - 2023-09-05

### 1.12.5 Changed

- Update AlphaFold to version v2.3.2

---

## [1.12.4] - 2023-08-29

### 1.12.4 Changed

- Removed VPC configuration for Lambda-based custom resources to better support bring-you-own-VPC

---

## [1.12.3] - 2023-08-8

### 1.12.3 Fixed

- Updated scipy version in JackHMMER container in response to Dependabot alert.

---

## [1.12.2] - 2023-08-7

### 1.12.2 Added

- Optional p4d compute environment and job queue
- Deployment script for updating existing stacks

### 1.12.2 Changed

- Updating the cfn stack will now restart the CodeBuild jobs for all module containers

---

## [1.12.1] - 2023-07-26

### 1.12.1 Fixed

- Fixed bug in miniconda version for OpenFold, ESMFold, and NextFlow container images.

### 1.12.1 Changed

- Updated ESMFold container to match Amazon HealthOmics implementation and improve run outputs.

---

## [1.12.0] - 2023-07-19

### 1.12.0 Fixed

- Due to a lack of ARM in the HMMER package, reverting all non-accelerated jobs back to x86 architecture.

---

## [1.11.0] - 2023-07-06

### 1.11.0 Added

- Support for running AlphaFold jobs on non-accelerated instance types by setting `gpu=0` and submitting job to the CPUOnDemandJobQueue.
- Support for passing a dictionary of environment variables to AlphaFold jobs.
- Support for c/m/r6i instance types in the CPU on-demand compute environment.
- Support for c/m/r7g instance types on the Graviton on-demand and spot compute environments.
- Added model version information to README.

---

## [1.10.4] - 2023-06-12

### 1.10.4 Changed

- Updated AlphaFold CUDA version to 11.6

### 1.10.4 Fixed

- Fixed bug causing hhsearch errors on g5 instance types

---

## [1.10.3] - 2023-06-07

### 1.10.3 Fixed

- Fixed bug preventing load of nested stacks.

---

## [1.10.2] - 2023-06-06

### 1.10.2 Changed

- Added additional error handling to batchfold package.

---

## [1.10.1] - 2023-05-31

### 1.10.1 Fixed

- Added paginators to boto3 calls in batchfold environment class to address issue with deploying multiple stacks.

---

## [1.10.0] - 2023-05-25

### 1.10.0 Added

- Added support for NextFlow

---

## [1.9.0] - 2023-04-21

### 1.9.0 Added

- Added support for RFDiffusion

### 1.9.0 Changed

- Removed deprecated notebook for RFDesign
- Update ProteinMPNN container and job definition to use accelerated instance types.
- Add new option to ProteinMPNN inference script to remove the input sequence from the output

### 1.9.0 Fixed

- Fixed bug with ESMFold metric logging for multiple sequence inputs.

---

## [1.8.0] - 2023-04-10

### 1.8.0 Added

- Added support for DiffDock
- Added ECR URIs as outputs to the CloudFormation template

### 1.8.0 Fixed

- Fixed typo in the pip install cell of the AlphaFold-Multimer notebook
- Updated version of py3Dmol to address structure visualization issue in JupyterLab

---

## [1.7.2] - 2023-03-21

### 1.7.2 Fixed

- Pin base image for ProteinMPNN and download containers.

---

## [1.7.1] - 2023-02-21

### 1.7.1 Fixed

- Check for AWSServiceRoleForEC2Spot and AWSServiceRoleForEC2SpotFleet service-linked roles before creating duplicates

---

## [1.7.0] - 2023-02-21

### 1.7.0 Added

- Updated AlphaFold version to 2.3.1
- Added support for OmegaFold release 2 models
- Moved reference data to S3 for faster stack provisioning

---

## [1.6.2] - 2023-02-10

### 1.6.2 Fixed

- Pinned commit of OmegaFold repository

---

## [1.6.1] - 2022-12-27

### 1.6.1 Fixed

- Fixed issue with ESMFold Dockerfile

---

## [1.6.0] - 2022-12-1

### 1.6.0 Added

- Added support for ProteinMPNN jobs.

### 1.6.0 Changed

- Incorporated ProteinMPNN into protein design notebook.
- Added list of supported algorithms to README index.

---

## [1.5.1] - 2022-11-30

### 1.5.1 Fixed

- Fixed ESMFold issue with malformed fasta headers
- Pinned ESMFold dependencies

---

## [1.5.0] - 2022-11-18

### 1.5.0 Added

- Added support for ESMFold jobs.

---

## [1.4.0] - 2022-11-18

### 1.4.0 Added

- Added support for RFDesign hallucination and inpainting jobs.

---

## [1.3.0] - 2022-10-12

### 1.3.0 Added

- Added default support for multi-AZ deployments.
- Added an S3 data repository to FSx for Lustre deployment to back up reference data.

### 1.3.0 Fixed

- Added exception handling for invalid JackHMMER `db_preset` values.

### 1.3.0 Changed

- Improved clarity of `prep_databases.py` script.
- Removed Batch compute environments for amd instance types.
- Moved S3 resource to root CloudFormation stack and updated retention policy
- Updated cost estimates.

---

## [1.2.0] - 2022-10-12

### 1.2.0 Added

- Support for the OmegaFold LLM-based protein folding algorithm.

### 1.2.0 Fixed

- Resolved issue with plotting protein structures on SageMaker Notebook instances.

### 1.2.0 Changed

- Update AlphaFold version to v2.2.4.
- Added nbhelper functions to BatchFold package
- Improved README

---

## [1.1.0] - 2022-09-19

### 1.1.0 Changed

- Refactor code to make algorithms more modular. This will make it easier to add new ones!

---

## [1.0.0] - 2022-09-08

### 1.0.0 Added

- Initial Release
