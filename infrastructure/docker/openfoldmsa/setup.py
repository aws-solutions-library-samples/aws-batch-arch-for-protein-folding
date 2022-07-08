from setuptools import setup, find_packages

setup(
    name='openfoldmsa',
    version='1.0.0',
    description='Code for generating multiple sequence alignment (MSA) files using the OpenFold implementation of Jackhmmer.',
    author='Brian Loyal',
    author_email='bloyal@amazon.com',
    license='Apache License, Version 2.0',
    url='https://github.com/aws-samples/aws-batch-architecture-for-alphafold',
    packages=find_packages(),
    python_requires='>=3',
    install_requires=[
        'biopython',
        'boto3',
        'datetime'
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