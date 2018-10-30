Save and load data from disk
----------------------------

AbPyTools uses its own custom object serialisation to save objects.

- FASTA

The simplest and most used format to save primary sequence data is the FASTA format.
Because it is so simple, it easy to parse and work with.
Here is an example of a FASTA file: ::
    >seq0
    FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF
    >seq1
    KYRTWEEFTRAAEKLYQADPMKVRVVLKYRHCDGNLCIKVTDDVVCLLYRTDQAQDVKKIEKFHSQLMRLME LKVTDNKECLKFKTDQAQEAKKMEKLNNIFFTLM

The sequence identifier is preceded by a ">" and the following line has the sequence.
This could be any sequence, such as DNA or protein, and is usually easy to distinguish due
to the different alphabets used to represent these.

The information that can be stored in a FASTA file is very limited, and
more sophisticated representations, such as sequence numbering of antibodies
is not possible.

- JSON
JSON is a simple format that resembles a Python dictionary, and
is easy to parse and read. However, it is inefficient in terms of disk usage.
AbPyTools uses this format to serialise objects, and these are human readable. ::
    {
        "ordered_names": ["test"],
        "test": {
            "numbering": [
                H1,
                H2,
                ...
                ]
            "numbering_scheme": "chothia",
            ...
    }

- Protocol Buffers
Protocol Buffers are relatively straight forward to use, and these
decode and encode messages in binary format, with a focus on speed
of encoding/decoding and disk usage reduction. This format is about
5 times more space efficient and at least as fast to read and write,
compared to JSON. The downside is that it is not human readable and
requires the precompiled parser to perform the encoding and decoding.
The message structure is defined in .proto files which can be found
in ``abpytools/core/formats``.