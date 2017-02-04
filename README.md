# AbPyTools
AbPyTools is a Python package to extract information from heavy and light antibody chain sequences. Using the built-in 
Antibody and AntibodyCollection it is very easy to manipulate the data and do more specific analysis with custom scripts.

This package is still in its early days and and lacks detailed documentation...!
Below are some existing features and some planned future additions.
At the moment it is updated several times a day until I am happy with the base functions. In the future the development 
will be performed in a separate branch.

AbPyTools features:
- 
- obtain Antibody numbering by querying AbNum (http://www.bioinf.org.uk/abs/abnum/)
- calculate hydrophobicity matrix
- calculate molecular weight
- calculate pI (isoelectric point)
- calculate CDR length
- higher level class that can load data from several antibody sequences
  - reads in data from FASTA format
  - calculates hydrophobicity matrix for whole dataset
  - get all the data already mentioned above
  - access CDR and framework sequences
- high level function to easily plot CDR length using a FASTA file as input
- currently only has Python 3 compatibility
- load and write antibody sequences in FASTA and json formats

Stuff that will be added/ worked on next:
- 
- Add remaining antibody numbering schemes
- Currently working on log-odds scores
- write tutorials
- adding some useful functions, such as comparing sequence with available datasets

- save Ab in new format to store all the stats that are calculated (i.e. numbering, MW, pI, ...). Most likely will store 
these in a file with .json format
- work on Python 2 and 3 compatibility
  - the main issue is the use of urllib (Python 3), instead of urllib2 (Python 2)
  - could use requests package instead 
- write tests
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