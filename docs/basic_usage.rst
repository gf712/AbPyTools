Basic usage
-----------

If this is your first time using AbPyTools you won't have any
JSON or pb2 files that are compatible with it's parsers.
The simplest is to either instantiate a :class:`Chain`
from a string, or load from a FASTA file with
:class:`ChainCollection`.

- :class:`Chain`


This object represents either a light or heavy chain of an
antibody. ::
    >>> from abpytools.core import Chain
    >>> my_seq = "MYFANCYLIGHTCHAIN"
    >>> chain = Chain.load_from_string(sequence=my_seq)

- :class:`ChainCollection`


This is a container that accesses :class:`Chain` objects and
used to perform analysis of a dataset. the easiest is to load
sequences from a string. ::
    >>> from abpytools.core import ChainCollection
    >>> collection = ChainCollection.load_from_fasta(path="secret.fasta")


Containers such as :class:`ChainCollection` can then be dumped
on disk to be used later on. This can be done either using
JSON format or pb2 (Protocol Buffer). ::
    >>> from abpytools.core import ChainCollection
    >>> collection = ChainCollection.load_from_fasta(path="secret.fasta")
    >>> # save to JSON
    >>> collection.save_to_json("loaded_secret")
    >>> # save to pb2
    >>> collection.save_to_pb2("loaded_secret")


And in a different session these two can be loaded just as
easily. ::
    >>> from abpytools.core import ChainCollection
    >>> collection_json = ChainCollection.load_from_json("loaded_secret.json")
    >>> collection_pb2 = ChainCollection.load_from_pb2("loaded_secret.pb2")