# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from setuptools import setup, find_packages

setup(
    name='batchfold',
    version='1.12.7',
    description='A modular architecture for running protein structure analysis on AWS Batch.',
    author='Brian Loyal',
    author_email='bloyal@amazon.com',
    license='Apache License, Version 2.0',
    url='https://github.com/aws-samples/aws-batch-architecture-for-alphafold',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3',
    install_requires=[
        'biopython',
        'boto3',
        'datetime',
        'requests'
    ],
    tests_require=[
        'pytest',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        "Programming Language :: Python :: 3",
    ],
)