from . import ChainProto


def get_protobuf_numbering_scheme(numbering_scheme):
    """
    Returns ChainProto field value of a given numbering scheme
    Args:
        numbering_scheme:

    Returns:

    """
    if numbering_scheme == "kabat":
        proto_numbering_scheme = ChainProto.KABAT
    elif numbering_scheme == "chothia":
        proto_numbering_scheme = ChainProto.CHOTHIA
    elif numbering_scheme == "chothia_ext":
        proto_numbering_scheme = ChainProto.CHOTHIA_EXT
    else:
        raise ValueError(f"{numbering_scheme} numbering scheme not supported by protobuf "
                         f"serialisation yet!")

    return proto_numbering_scheme


def get_numbering_scheme_from_protobuf(proto_numbering_scheme):
    """
    Returns ChainProto field value of a given numbering scheme
    Args:
        numbering_scheme:

    Returns:

    """
    if proto_numbering_scheme == ChainProto.KABAT:
        numbering_scheme = "kabat"
    elif proto_numbering_scheme == ChainProto.CHOTHIA:
        numbering_scheme = "chothia"
    elif proto_numbering_scheme == ChainProto.CHOTHIA_EXT:
        numbering_scheme = "chothia_ext"
    else:
        raise ValueError(f"ChainProto numbering scheme {proto_numbering_scheme} is not compatible with abpytools!")

    return numbering_scheme
