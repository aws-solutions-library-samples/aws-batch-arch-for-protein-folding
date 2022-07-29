# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2022-07-29

### Added
- Added support for MMSeqs2 and OpenFold algorithms
- New BatchFold library for easier job submission

### Changed
- Significant refactoring of code and tests

## [1.0.4] - 2022-06-24

### Changed
- Added support for AlphaFold v2.2.2


## [1.0.3] - 2022-06-03

### Fixed
- Resolved issue where timings.json file was not properly passed between the data prep and folding jobs (Identified by PR https://github.com/aws-samples/aws-batch-architecture-for-alphafold/pull/3)
- Resolved compatibility issue with Protobuf version (https://github.com/deepmind/alphafold/issues/478)

### Changed

- Added support for AlphaFold v2.2.0
- Added support for existing VPCs and FSx for Lustre instances
- Added option to automatically download reference data to FSx
- Decreased data prep job cost by as much as 75% by using spot instance types.
- Updated data download script to pull sequence databases and parameters from S3
- Decreased data download time by 41% by retrieving ref data in parallel
- Refactored CloudFormation template to use nested stacks
- Added CloudFormation tests with Taskcat
- Updated folding container to Python 3.9, Ubuntu 20.04, and mamba

## [1.0.2] - 2022-05-20

### Fixed

### Changed

- Added additional download scripts to help address the instability of the public PDB mmcif mirror.
- Moved download jobs into private subnet for increased security.

## [1.0.1] - 2022-02-25

### Fixed

- Fix bug in download-Dockerfile that was causing download jobs to run longer than needed.
- Fix typo in section headers in AWS-AlphaFold notebook.

### Changed

- Add table of contents to AWS-AlphaFold notebook.

## [1.0] - 2022-02-24

### Changed

- Open Source Release

## [0.2] - 2022-02-23

### Changed

- Documentation updates
- Security and IP updates

## [0.1] - 2022-02-10

### Added

- Initial Release
- Submit 1- and 2-step protein-folding jobs to AWS Batch using a Jupyter Notebook.