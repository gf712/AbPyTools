from abpytools.core.flags import numbering, pi, backend, options, hydrophobicity, chain, formats

# create alias to flag modules
NUMBERING_FLAGS = numbering
CHAIN_FLAGS = chain
OPTION_FLAGS = options
PI_FLAGS = pi
BACKEND_FLAGS = backend
HYDROPHOBICITY_FLAGS = hydrophobicity
FORMAT_FLAGS = formats

# from Python 3.7 can do this
# import abpytools.core.flags.numbering as NUMBERING_FLAGS
# import abpytools.core.flags.chain as CHAIN_FLAGS
# import abpytools.core.flags.options as OPTION_FLAGS
# import abpytools.core.flags.pi as PI_FLAGS
# import abpytools.core.flags.hydrophobicity as HYDROPHOBICITY_FLAGS
# import abpytools.core.flags.backend as BACKEND_FLAGS

__all__ = ["NUMBERING_FLAGS",
           "CHAIN_FLAGS",
           "OPTION_FLAGS",
           "PI_FLAGS",
           "HYDROPHOBICITY_FLAGS",
           "BACKEND_FLAGS",
           "FORMAT_FLAGS"]
