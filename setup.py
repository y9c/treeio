#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open("README.md") as readme_file:
    README = readme_file.read()

REQUIREMENTS = [
    "Click>=6.0",
]

SETUP_REQUIREMENTS = [
    "pytest-runner",
]

TEST_REQUIREMENTS = [
    "pytest",
]

setup(
    author="Ye Chang",
    author_email="yech1990@gmail.com",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="A pure python parser for tree data format IO",
    entry_points={"console_scripts": ["treeio=treeio.cli:cli",],},
    install_requires=REQUIREMENTS,
    license="MIT license",
    long_description=README,
    include_package_data=True,
    keywords="treeio",
    name="treeio",
    packages=find_packages(include=["treeio"]),
    setup_requires=SETUP_REQUIREMENTS,
    test_suite="tests",
    tests_require=TEST_REQUIREMENTS,
    url="https://github.com/yech1990/treeio",
    version="0.0.0.dev0",
    zip_safe=False,
)
