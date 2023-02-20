# Original Copyright 2021 DeepMind Technologies Limited
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Install script for setuptools."""

from setuptools import find_packages
from setuptools import setup

setup(
    name='alphafold-jackhmmer',
    version='1.0',
    description='An implementation of the Jackhmmer MSA creation pipeline, based off the data pipeline included in AlphaFold v2.3.0',
    author='AWS',
    license='Apache License, Version 2.0',
    packages=find_packages(),
    install_requires=[
        'absl-py',
        'biopython',
        'numpy',
        'pandas',
        'scipy',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
)
