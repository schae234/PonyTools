#!/usr/bin/env python3

import os
import io
import re

from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools import setup, find_packages, Extension

from subprocess import check_call

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

    description = 'Computational tools that you will ride into the FUTURE.',
    url = 'https://github.com/schae234/PonyTools',

    author = 'Rob Schaefer',
    author_email = 'rob@linkage.io',
    license = "Copyright Linkage Analytics 2016. Available under the MIT License",

    classifiers=[
	# How mature is this project? Common values are
	#   3 - Alpha
	#   4 - Beta
	#   5 - Production/Stable
	'Development Status :: 4 - Beta',

	# Indicate who your project is intended for
	'Intended Audience :: Developers',
	'Topic :: Software Development :: Build Tools',

	# Pick your license as you wish (should match "license" above)
	 'License :: OSI Approved :: MIT License',

	# Specify the Python versions you support here. In particular, ensure
	# that you indicate whether you support Python 2, Python 3 or both.
	'Programming Language :: Python :: 3',
	'Programming Language :: Python :: 3.6',
    ],
    keywords='data genetics biology', 
    project_urls={
        'Source' : 'https://github.com/schae234/PonyTools',
        'Tracker' : 'https://github.com/schae234/PonyTools/issues'
    },

    packages = find_packages(),
    scripts = [
       'ponytools/cli/ponytools'     
    ],
    include_package_data=True,
    ext_modules = [],
    cmdclass = {
        #'build_ext': build_ext
    },
    package_data={
        'MNEc2MAnnot':'data/MNEc2M_Annotation.csv.gz', 
        'ImputationMakefile':'scripts/ImputationMakefile'
    },
    dependency_links = [
        ''
    ],
    install_requires = [
        'matplotlib>=1.4.3',
        'pandas>=0.16.2',
        'ipython>=4.1.2',
        'numpy>=1.12.0',
        'scipy>=0.16.0',
        'ipdb>=0.8',
        'cython>=0.25.2',
        'tqdm>=4.23.2',
        'pysam>=0.14.1',
        'minus80>=0.1.2',
        'locuspocus>=0.1.1'
    ],


)
