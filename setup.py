import sys
from setuptools import setup, find_packages
from setuptools import Extension
import os

SETUP_METADATA = \
               {
    "name": "phamb",
    "description": "Phages from metagenomic binning",
    "url": "https://github.com/RasmussenLab/phamb",
    "version": "1.0.1",
    "license": "MIT",
    "packages": ['phamb'],
    "package_data": {'phamb': ['dbs/RF_model.python39.sav']},
    "python_requires": ">=3.9",
    "install_requires": ["scikit-learn=1.0.2"],
    "classifiers":[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    "scripts":['phamb/run_RF.py','phamb/split_contigs.py','phamb/vambtools.py','phamb/run_RF_modules.py']
    }

setup(**SETUP_METADATA)