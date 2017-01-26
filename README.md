# AbPyTools
Repo with some scripts to analyse antibody sequences

This repo is still in its early days so some scripts will not be very user friendly...!
Over time I will rewrite these to be easier to use and add some data (i.e. numbering schemes for Ab)

What it can do now:
- 
- obtain Antibody numbering by querying AbNum
- calculate hydrophobicity matrix
- calculate molecular weight
- calculate pI (isoelectric point)
- calculate CDR length

Stuff that will be added next:
- 
- class that holds multiple Ab objects
- adding some useful functions, such as comparing sequence with available datasets
- Ab loader (from FASTA format probably)
- write out Ab sequences (FASTA format)
- save Ab in new format to store all the stats that are calculated (i.e. numbering, MW, pI, ...). Most likely will store 
these in a file with .json format
- work on Python 2 and 3 compatibility
- write setup.py
- write tests
