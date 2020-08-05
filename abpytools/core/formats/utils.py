from ..flags import *
if BACKEND_FLAGS.HAS_PROTO:
    from . import ChainProto

def get_protobuf_numbering_scheme(numbering_scheme):
    """
    Returns ChainProto field value of a given numbering scheme
    Args:
        numbering_scheme:

    Returns:

    """
    if numbering_scheme == NUMBERING_FLAGS.KABAT:
        proto_numbering_scheme = ChainProto.KABAT
    elif numbering_scheme == NUMBERING_FLAGS.CHOTHIA:
        proto_numbering_scheme = ChainProto.CHOTHIA
    elif numbering_scheme == NUMBERING_FLAGS.CHOTHIA_EXT or numbering_scheme == NUMBERING_FLAGS.MARTIN:
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
        numbering_scheme = NUMBERING_FLAGS.KABAT
    elif proto_numbering_scheme == ChainProto.CHOTHIA:
        numbering_scheme = NUMBERING_FLAGS.CHOTHIA
    elif proto_numbering_scheme == ChainProto.CHOTHIA_EXT:
        numbering_scheme = NUMBERING_FLAGS.CHOTHIA_EXT
    else:
        raise ValueError(f"ChainProto numbering scheme {proto_numbering_scheme} is not compatible with abpytools!")

    return numbering_scheme
