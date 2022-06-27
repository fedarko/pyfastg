# NOTE: Derived from https://github.com/biocore/qurro/blob/master/setup.py

from setuptools import find_packages, setup

classes = """
    Development Status :: 3 - Alpha
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
"""
classifiers = [s.strip() for s in classes.split("\n") if s]

description = "Minimal Python library for parsing SPAdes FASTG files"

with open("README.md") as f:
    long_description = f.read()

setup(
    name="pyfastg",
    version="0.1.0",
    license="MIT",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Nick Waters, Marcus Fedarko",
    maintainer="Marcus Fedarko",
    maintainer_email="mfedarko@ucsd.edu",
    url="https://github.com/fedarko/pyfastg",
    classifiers=classifiers,
    packages=find_packages(),
    install_requires=["networkx >= 2"],
    extras_require={"dev": ["pytest", "pytest-cov", "flake8", "black"]},
)
