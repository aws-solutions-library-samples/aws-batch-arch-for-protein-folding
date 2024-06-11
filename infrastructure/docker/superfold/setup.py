from setuptools import setup, find_packages

setup(name='superfold',
      version='1.0.0',
      description='SuperFold: Fast single-sequence predictor for IPD labs.',
      author='UW Baker Lab',
      url='https://github.com/rdkibler/superfold.git',
      scripts=["scripts/run_superfold.py"],
      packages=find_packages(),
      install_requires=['pytorch'])