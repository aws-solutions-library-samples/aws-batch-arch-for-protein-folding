# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.9.0] - 2023-04-21

### 1.9.0 Added

- Added support for RFDiffusion

### 1.9.0 Changed

- Removed deprecated notebook for RFDesign
- Update ProteinMPNN container and job definition to use accelerated instance types.
- Add new option to ProteinMPNN inference script to remove the input sequence from the output

### 1.8.0 Fixed

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
