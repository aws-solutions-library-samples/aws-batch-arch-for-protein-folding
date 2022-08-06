# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
from setuptools import setup, find_packages

from torch.utils.cpp_extension import BuildExtension, CUDAExtension

version_dependent_macros = [
    "-DVERSION_GE_1_1",
    "-DVERSION_GE_1_3",
    "-DVERSION_GE_1_5",
]

extra_cuda_flags = [
    "-std=c++14",
    "-maxrregcount=50",
    "-U__CUDA_NO_HALF_OPERATORS__",
    "-U__CUDA_NO_HALF_CONVERSIONS__",
    "--expt-relaxed-constexpr",
    "--expt-extended-lambda",
]

compute_capabilities = set(
    [
        (7, 0),  # V100 for p3
        (7, 5),  # T4 for g4dn
        (8, 0),  # A100 for p4dn
        (8, 6),  # A10 for g5
    ]
)

cc_flag = []
for major, minor in list(compute_capabilities):
    cc_flag.extend(
        [
            "-gencode",
            "arch=compute_{major}{minor},code=sm_{major}{minor}".format(major=major, minor=minor),
        ]
    )

extra_cuda_flags += cc_flag

setup(
    name="openfold",
    version="1.0.0",
    description="A PyTorch reimplementation of DeepMind's AlphaFold 2",
    author="Gustaf Ahdritz & DeepMind",
    author_email="gahdritz@gmail.com",
    license="Apache License, Version 2.0",
    url="https://github.com/aqlaboratory/openfold",
    packages=find_packages(exclude=["tests", "scripts"]),
    include_package_data=True,
    package_data={
        "openfold": ["utils/kernel/csrc/*"],
        "": ["resources/stereo_chemical_props.txt"],
    },
    ext_modules=[
        CUDAExtension(
            name="attn_core_inplace_cuda",
            sources=[
                "openfold/utils/kernel/csrc/softmax_cuda.cpp",
                "openfold/utils/kernel/csrc/softmax_cuda_kernel.cu",
            ],
            include_dirs=[
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "openfold/utils/kernel/csrc/",
                )
            ],
            extra_compile_args={
                "cxx": ["-O3"] + version_dependent_macros,
                "nvcc": (
                    ["-O3", "--use_fast_math"]
                    + version_dependent_macros
                    + extra_cuda_flags
                ),
            },
        )
    ],
    cmdclass={"build_ext": BuildExtension},
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
)
