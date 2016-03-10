#!/usr/bin/env python3

import os
import io
import re

from setuptools import setup, find_packages, Extension

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name = 'ponytools',
    version = find_version('ponytools','__init__.py'),
    packages = find_packages(),
    scripts = [
       'ponytools/cli/ponytools'     
    ],

    install_requires = [
        'matplotlib>=1.4.3',
        'pandas>=0.16',
        'ipython>=4.1.2',
        'scipy>=0.17.0',
        'ipdb>=0.8'
    ],

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Computational tools that you will ride into the FUTURE.',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/PonyTools',

)
