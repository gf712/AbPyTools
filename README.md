# AbPyTools
Repo with some scripts to analyse antibody sequences

This repo is still in its early days so some scripts will not be very user friendly...!
Over time I will rewrite these to be easier to use and add some data (i.e. numbering schemes for Ab)

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
- currently only has Python 3 compatibility

Stuff that will be added/ worked on next:
- 
- write tutorials
- adding some useful functions, such as comparing sequence with available datasets
- write out Ab sequences (e.g. FASTA format)
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