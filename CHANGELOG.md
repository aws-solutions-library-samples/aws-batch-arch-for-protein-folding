# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
