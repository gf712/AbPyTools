from abpytools.formats.utils import get_protobuf_numbering_scheme, get_numbering_scheme_from_protobuf
from abpytools.core.chain import Chain


def _name_wrapper(func):
    """
    Internal usage only counter decorator to format sequence names.
    Args:
        func:

    Returns:

    """
    def helper(*args):
        helper.idi += 1
        name = func(*args)
        if isinstance(name, str):
            name = name.replace("%%IDI%%", str(helper.idi))
        else:
            names = list()
            for name_i in name:
                names.append(name_i.replace("%%IDI%%", str(helper.idi)))
            name = names
        return name
    helper.idi = 0
    return helper


@_name_wrapper
def _chain_name_getter(antibody_object):
    """
    Internal function that sets the name of Chain objects.
    This is required because it is possible to work with a ChainCollection
    where a Chain does not have a name. However, a name is
    required internally to store data on disk.

    Args:
        antibody_object (Chain):

    Returns:

    """
    # if antibody does not have name, generate name:
    # ID_chain_idi, where chain is heavy/light, idi is i = [1,..,N]
    if len(antibody_object.name) > 0:
        name = antibody_object.name
    else:
        name = f"ID_{antibody_object.chain}_%%IDI%%"
    return name


@_name_wrapper
def _fab_name_getter(name_i, light_chain, heavy_chain):
    """
    Internal function that sets the name of Fab objects
    when they are stored on disk.
    Args:
        name_i (str):
        light_chain (Chain):
        heavy_chain (Chain):

    Returns:

    """
    # if antibody does not have name, generate name:
    # ID_chain_idi, where chain is heavy/light, idi is i = [1,..,N]
    if len(name_i) > 0:
        names = name_i
    else:
        names = []
        for chain in [light_chain, heavy_chain]:
            if len(chain.name) > 0:
                names.append(chain.name)
            else:
                names.append(f"ID_{chain.chain}_%%IDI%%")
    return names


##############################################################################
#                                   SAVERS
##############################################################################

def json_FabCollection_formatter(fab_object):
    """
    Internal function to serialise FabCollection objects in JSON format
    Args:
        fab_object (FabCollection):

    Returns:

    """

    light_antibody_json = json_ChainCollection_formatter(fab_object._light_chains)
    heavy_antibody_json = json_ChainCollection_formatter(fab_object._heavy_chains)

    fab_data = dict()
    ordered_names = []

    for light_antibody, heavy_antibody in zip(
            light_antibody_json['ordered_names'], heavy_antibody_json['ordered_names']):
        name = f"{light_antibody_json[light_antibody]['name']}-" \
               f"{heavy_antibody_json[heavy_antibody]['name']}"
        fab_data[name] = [light_antibody_json[light_antibody],
                          heavy_antibody_json[heavy_antibody]]
        ordered_names.append(name)
    fab_data['ordered_names'] = ordered_names

    return fab_data


def json_ChainCollection_formatter(chain_objects):
    """
    Internal function to serialise ChainCollection objects in JSON format
    Args:
        chain_objects (ChainCollection):

    Returns:

    """
    data = dict()
    ordered_names = []
    # reset wrapper counter
    _chain_name_getter.idi = 0
    for antibody in chain_objects:
        name = _chain_name_getter(antibody)
        ordered_names.append(name)
        antibody_dict = antibody.ab_format()
        # use name generated from _chain_name_getter instead
        antibody_dict['name'] = name
        data[name] = antibody_dict
    data['ordered_names'] = ordered_names
    return data


def pb2_ChainCollection_formatter(chain_objects, proto_parser, reset_status=True):
    """
    Internal function to serialise a ChainCollection object to .pb2 format
    according to definitition in 'format/chain.proto'.

    Args:
        chain_objects (ChainCollection):
        proto_parser (ChainCollectionProto):
        reset_status (bool):

    Returns:

    """
    # reset wrapper counter
    if reset_status:
        _chain_name_getter.idi = 0
    for chain_object in chain_objects:
        pb2_add_chain(chain_object, proto_parser)


def pb2_FabCollection_formatter(fab_object, proto_parser, reset_status=True):
    """
    Internal function to serialise a FabCollection object to .pb2 format
    according to definitition in 'format/fab.proto'.

    Args:
        fab_object (FabCollection):
        proto_parser (FabCollectionProto):
        reset_status (bool):

    Returns:

    """
    # reset wrapper counter
    if reset_status:
        _chain_name_getter.idi = 0
    for name_i, light_chain_i, heavy_chain_i in zip(fab_object.names,
                                                    fab_object._light_chains,
                                                    fab_object._heavy_chains):

        proto_fab = proto_parser.fabs.add()
        name_i = _fab_name_getter(name_i, light_chain_i, heavy_chain_i)
        proto_fab.name = name_i
        add_Chain_to_protobuf(light_chain_i, proto_fab.light_chain)
        add_Chain_to_protobuf(heavy_chain_i, proto_fab.heavy_chain)


def pb2_add_chain(chain_object, proto_parser):
    """
    Populates a protobuf ProtoChain message from Chain object
    and adds it to ChainCollectionProto
    Args:
        chain_object (Chain):
        proto_parser (ChainCollectionProto):

    Returns:

    """
    proto_antibody = proto_parser.chains.add()
    add_Chain_to_protobuf(chain_object, proto_antibody)


def add_Chain_to_protobuf(antibody_obj, proto_obj):
    """
    Helper function to populate a ProtoChain message.
    Args:
        antibody_obj (Chain):
        proto_obj (ChainProto):

    Returns:

    """
    proto_obj.name = _chain_name_getter(antibody_obj)
    proto_obj.sequence = antibody_obj.sequence
    proto_obj.numbering_scheme = get_protobuf_numbering_scheme(antibody_obj.numbering_scheme)
    proto_obj.numbering.extend(antibody_obj.numbering)
    proto_obj.chain_type = antibody_obj.chain
    proto_obj.MW = antibody_obj.mw
    for x in {**antibody_obj.ab_regions()[0], **antibody_obj.ab_regions()[1]}.items():
        proto_obj_cdr = proto_obj.region.add()
        proto_obj_cdr.region_name = x[0]
        proto_obj_cdr.region_positions.extend(x[1])
    proto_obj.pI = antibody_obj.pI


##############################################################################
#                                   PARSERS
##############################################################################
def fasta_ChainCollection_parser(raw_fasta, numbering_scheme):

    names = list()
    sequences = list()
    antibody_objects = list()
    for line in raw_fasta:
        if line.startswith(">"):
            names.append(line.replace("\n", "")[1:])
        # if line is empty skip line
        elif line.isspace():
            pass
        else:
            sequences.append(line.replace("\n", ""))
    if len(names) != len(sequences):
        raise ValueError("Error reading file: make sure it is FASTA format")

    for name, sequence in zip(names, sequences):
        antibody_objects.append(Chain(name=name, sequence=sequence,
                                      numbering_scheme=numbering_scheme))

    return antibody_objects


def pb2_FabCollection_parser(proto_parser):
    from abpytools.core.fab import Fab
    fab_objects = list()
    for fab_i in proto_parser.fabs:

        light_chain = pb2_Chain_parser(fab_i.light_chain)
        heavy_chain = pb2_Chain_parser(fab_i.heavy_chain)

        fab = Fab(light_chain=light_chain, heavy_chain=heavy_chain, name=fab_i.name)

        fab_objects.append(fab)

    return fab_objects


def pb2_ChainCollection_parser(proto_parser):
    antibody_objects = list()
    names = list()
    for chain_i in proto_parser.chains:

        antibody_i = pb2_Chain_parser(chain_i)

        antibody_i._loading_status = 'Loaded'
        names.append(chain_i.name)
        antibody_objects.append(antibody_i)

    return antibody_objects


def pb2_Chain_parser(proto_chain):
    """
    Populate Chain object from protobuf file

    Args:
        proto_chain (ChainProto):

    Returns:

    """
    chain_obj = Chain(name=proto_chain.name, sequence=proto_chain.sequence)

    if proto_chain.numbering_scheme:
        chain_obj._numbering_scheme = get_numbering_scheme_from_protobuf(proto_chain.numbering_scheme)

    if proto_chain.numbering:
        # convert google.protobuf.pyext._message.RepeatedScalarContainer to list
        chain_obj.numbering = list(proto_chain.numbering)
    else:
        chain_obj.numbering = chain_obj.ab_numbering()

    if proto_chain.chain_type:
        chain_obj._chain = proto_chain.chain_type
    else:
        chain_obj._chain = Chain.determine_chain_type(proto_chain.numbering)

    if proto_chain.MW:
        chain_obj.mw = proto_chain.MW
    else:
        chain_obj.mw = chain_obj.ab_molecular_weight()

    if proto_chain.region:
        regions = dict()
        for region in proto_chain.region:
            regions[region.region_name] = region.region_positions
        chain_obj.CDR = regions
    else:
        chain_obj.CDR = chain_obj.ab_regions()

    if proto_chain.pI:
        chain_obj.pI = proto_chain.pI
    else:
        chain_obj.pI = chain_obj.ab_pi()

    return chain_obj


def json_FabCollection_parser(raw_data):
    from abpytools.core.fab import Fab
    fab_objects = list()

    ordered_names = raw_data.pop('ordered_names')

    for key_i in ordered_names:
        light_chain_i = json_Chain_parser(raw_data[key_i][0], raw_data[key_i][0]["name"])
        heavy_chain_i = json_Chain_parser(raw_data[key_i][1], raw_data[key_i][1]["name"])
        fab_i = Fab(light_chain=light_chain_i, heavy_chain=heavy_chain_i, name=key_i)
        fab_objects.append(fab_i)

    return fab_objects


def json_ChainCollection_parser(raw_data):

    antibody_objects = list()

    ordered_names = raw_data.pop('ordered_names')

    for key_i in ordered_names:
        chain_i = json_Chain_parser(raw_data[key_i], key_i)
        antibody_objects.append(chain_i)

    return antibody_objects


def json_Chain_parser(antibody_dict, name):
    antibody_i = Chain(name=name, sequence=antibody_dict['sequence'])
    antibody_i.numbering = antibody_dict['numbering']
    antibody_i._chain = antibody_dict['chain']
    antibody_i.mw = antibody_dict['MW']
    antibody_i.CDR = antibody_dict["CDR"]
    antibody_i._numbering_scheme = antibody_dict["numbering_scheme"]
    antibody_i.pI = antibody_dict["pI"]
    antibody_i._loading_status = 'Loaded'
    return antibody_i
