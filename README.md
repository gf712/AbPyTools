AbPyTools
=========

[![Coverage Status](https://coveralls.io/repos/github/gf712/AbPyTools/badge.svg?branch=master)](https://coveralls.io/github/gf712/AbPyTools?branch=master)
[![Build Status](https://travis-ci.org/gf712/AbPyTools.svg?branch=master)](https://travis-ci.org/gf712/AbPyTools)
[![Code Health](https://landscape.io/github/gf712/AbPyTools/master/landscape.svg?style=flat)](https://landscape.io/github/gf712/AbPyTools/master)
[![Documentation Status](https://readthedocs.org/projects/abpytools/badge/?version=latest)](https://abpytools.readthedocs.io/en/latest/?badge=latest)


AbPyTools is a Python 3 package to extract information from heavy and light antibody chain sequences. Using the built-in
Antibody and ChainCollection it is very easy to manipulate the data and do more specific analysis with custom scripts.

This package is still in its early days and and lacks detailed documentation...!
Below are some existing features and some planned future additions.
At the moment it is updated several times a day until I am happy with the base functions. In the future the development 
will be performed in a separate branch.

AbPyTools features
- 
- obtain Antibody numbering by querying AbNum (http://www.bioinf.org.uk/abs/abnum/)
- analysis of scFv sequences with optimised backend
- higher level class that can load data from several antibody sequences
  - load and write antibody sequences in FASTA and json formats
  - calculates hydrophobicity matrix for whole dataset
  - get all the data already mentioned above
  - access CDR and framework sequences
- Work with heavy and light chains or combinations 
- high level function to easily plot CDR length using a FASTA file as input

Stuff that will be added/ worked on next
- 
- Add remaining antibody numbering schemes
- write tutorials
- adding some useful functions, such as comparing sequence with available datasets
- write high level code for more specific analysis
  - plot CDR lengths of antibodies
  
Cython
-
From version 0.3 AbPyTools will start using Cython to speed up numerical manipulations.
In the front end nothing will change, but installation from source will require Cython!

This new feature will speed up most calculations significantly. The backend uses code written 
from scratch that mimics numpy behaviour but runs much faster, since it is more specialised and lightweight.

Installing abpytools
-

## From source

Clone code from the GitHub repository

`git clone https://github.com/gf712/AbPyTools.git`

Change to package directory in your local machine

`cd path/to/AbPyTools`

Install package

`python setup.py install`

Run tests (recommended)

`python setup.py test`

## From pypi

`pip install abpytools` 

Import to python
-
`import abpytools`

Changelog
-
### v0.3.2 (release date 31/10/2018)
 - Added docs page
 - Fixed Pypi build with Cython files
### v0.3.1 (release date 30/10/2018)
 - Major:
   - Protobuf support to serialise core objects (`ChainCollection` and `FabCollection`)
     - speed up in saving and loading large files
     - files are about 5 times smaller
   - Dropped support for Python 3.5 to start using f-strings
   - API changes:
     - Cleaned up object instantiation and added factory functions (**This will break some old code but provides a cleaner iterface**)
 - Minor:
   - added Python 3.7 to Travis script (also dropped Python 3.5)
### v0.3 (release date 11/10/2018):
 - Implementation of backend calculations with Cython leading to speedups of several orders of magnitude
