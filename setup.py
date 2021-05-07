# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 11:37:46 2019

@author: zhang
"""


import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="cfDNApipe",
    version="1.0.3",
    author="Wei Zhang, Jiaqi Huang, Shuying He, Juhong Liu, Yu Liu",
    author_email="w-zhang16@mail.tsinghua.edu.cn",
    description="An Intergrated Pipeline For cfDNA Sequencing Data Analysis",
    license="Please see LICENSE.txt.",
    keywords=["cell free DNA", "WGS", "WGBS", "Fragmentation", "Methylation", "Virus", "SNV", "CNV"],
    url="https://xwanglabthu.github.io/cfDNApipe/",
    packages=find_packages(),
    package_data={"cfDNApipe": ["data/*", "temp/*"]},
    long_description=read("README.rst"),
    platforms="Linux/Unix",
)
