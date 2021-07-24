# _*_ coding: utf-8 _*_
"""
Time:     2021/7/22 18:37
Author:   WANG Bingchen
Version:  V 0.1
File:     setup.py
Describe: 
"""

import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="scend",
    version="1.0.0",
    author="Bingchen Wang",
    author_email="wangbingchen@buaa.edu.cn",
    description="SCEnd: A matrix factorization method for single-cell",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/deepomicslab/SCEnd",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
