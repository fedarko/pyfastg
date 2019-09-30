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

version = "0.0.0"

setup(
    name="pyfastg",
    version=version,
    description=description,
    long_description=description,
    author="Nick Waters",
    url="https://github.com/nickp60/pyfastg",
    classifiers=classifiers,
    packages=find_packages(),
    install_requires=["networkx"],
    extras_require={"dev": ["pytest", "pytest-cov", "flake8", "black"]},
)
