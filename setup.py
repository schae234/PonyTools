#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

setup(
    name = 'ponytools',
    version = '0.1.7',
    packages = find_packages(),
    scripts = [
       'bin/pyDi'     
    ],

    install_requires = [
        'matplotlib>=1.4.3',
        'pandas>=0.16',
        'ipdb>=0.8'
    ],

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Computational tools that you will ride into the FUTURE.',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/PonyTools',


)
