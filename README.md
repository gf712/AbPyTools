AbPyTools
=========

[![Coverage Status](https://coveralls.io/repos/github/gf712/AbPyTools/badge.svg?branch=master)](https://coveralls.io/github/gf712/AbPyTools?branch=master)
[![Build Status](https://travis-ci.org/gf712/AbPyTools.svg?branch=master)](https://travis-ci.org/gf712/AbPyTools)
[![Code Health](https://landscape.io/github/gf712/AbPyTools/master/landscape.svg?style=flat)](https://landscape.io/github/gf712/AbPyTools/master)

AbPyTools is a Python 3 package to extract information from heavy and light antibody chain sequences. Using the built-in
Antibody and ChainCollection it is very easy to manipulate the data and do more specific analysis with custom scripts.

This package is still in its early days and and lacks detailed documentation...!
Below are some existing features and some planned future additions.
At the moment it is updated several times a day until I am happy with the base functions. In the future the development 
will be performed in a separate branch.

AbPyTools features:
- 
- obtain Antibody numbering by querying AbNum (http://www.bioinf.org.uk/abs/abnum/)
- calculate molecular weight, pI (isoelectric point) and extinction coefficients
- calculate hydrophobicity matrix and CDR length
- higher level class that can load data from several antibody sequences
  - load and write antibody sequences in FASTA and json formats
  - calculates hydrophobicity matrix for whole dataset
  - get all the data already mentioned above
  - access CDR and framework sequences
- Work with heavy and light chains or combinations 
- high level function to easily plot CDR length using a FASTA file as input

Stuff that will be added/ worked on next:
- 
- Add remaining antibody numbering schemes
- write tutorials
- adding some useful functions, such as comparing sequence with available datasets
- write high level code for more specific analysis
  - plot CDR lengths of antibodies

Installing abpytools:
-
Clone code from the GitHub repository

`git clone https://github.com/gf712/AbPyTools.git`

Change to package directory in your local machine

`cd path/to/AbPyTools`

Install package

`python setup.py install`

Import to python
-
`import abpytools`