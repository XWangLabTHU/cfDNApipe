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
    name = "cfDNApipe",
    version = "0.0.4",
    author = "Wei Zhang",
    author_email = "w-zhang16@mail.tsinghua.edu.cn",
    description = "An Intergrated Pipeline For cfDNA Sequencing Data",
    license = "Please see LICENSE.txt.",
    keywords = ['cfDNA', 'WGS', 'WGBS'],
    url = "https://github.com/Honchkrow/cfDNApipe",
    packages = find_packages(),
	package_data = {'cfDNApipe': ['data/*']},
    long_description = read('README.rst'),
    platforms = "Linux/Unix, macOS"
)