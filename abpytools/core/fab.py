from .chain_collection import Chain, ChainCollection
import numpy as np
from .chain import calculate_charge
from abpytools.utils import DataLoader
from .helper_functions import germline_identity_pd, to_numbering_table


class Fab:
    """
    
    """
    def __init__(self, heavy_chain=None, light_chain=None, load=True, name=None):
        # As a convention the order will always be light chain (index 0) and then heavy chain (index 1)

        # check if it's an Chain class
        if heavy_chain is None and light_chain is None:
            raise ValueError('heavy_chain and light_chain must be provided')

        if isinstance(heavy_chain, Chain):
            self._heavy_chain = heavy_chain

        else:
            raise ValueError('heavy_chain must be a Chain object, but got {} instead'.format(type(heavy_chain)))

        if isinstance(light_chain, Chain):
            self._light_chain = light_chain

        else:
            raise ValueError('light_chain must be a Chain object, but got {} instead'.format(type(light_chain)))

        if self._heavy_chain.chain != 'heavy':
            raise ValueError("heavy_chain is not a heavy chain, it is {}".format(self._heavy_chain.chain))
        if self._light_chain.chain != 'light':
            raise ValueError("light_chain is not a light chain, it is {}".format(self._light_chain.chain))

        self._pair_sequence = self[0].sequence + self[1].sequence

        if isinstance(name, str):
            self._name = name

        elif name is None:
            self._name = 'ID1'

        # keep the name of the heavy and light chains internally to keep everything in the right order
        self._internal_heavy_name = self[1].name
        self._internal_light_name = self[0].name

    def load(self):
        self._heavy_chain.load()
        self._light_chain.load()

    def molecular_weight(self, monoisotopic=False):

        return self[0].ab_molecular_weight(monoisotopic=monoisotopic) +\
               self[1].ab_molecular_weight(monoisotopic=monoisotopic)

    def extinction_coefficient(self, reduced=False, normalise=False, **kwargs):

        light_ec = self[0].ab_ec(reduced=reduced, **kwargs)
        heavy_ec = self[1].ab_ec(reduced=reduced, **kwargs)

        if normalise:
            return (heavy_ec + light_ec) / (self.molecular_weight(**kwargs))
        else:
            return heavy_ec + light_ec

    def hydrophobicity_matrix(self, **kwargs):

        return np.concatenate((self[0].ab_hydrophobicity_matrix(**kwargs),
                               self[1].ab_hydrophobicity_matrix(**kwargs)))

    def charge(self, **kwargs):

        return np.concatenate((self[0].ab_charge(**kwargs), self[1].ab_charge(**kwargs)))

    def total_charge(self, ph=7.4, pka_database='Wikipedia'):

        available_pi_databases = ["EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell", "Wikipedia", "Lehninger",
                                  "Grimsley"]
        assert pka_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pka_database,
                                                                                    ', '.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties', data=['pI', pka_database])
        pka_data = data_loader.get_data()

        return calculate_charge(sequence=self.sequence, ph=ph, pka_values=pka_data)

    # def load_igblast_query(self, file_path, chain):
    #
    #     if chain.lower() == 'light':
    #         self._light_chains.load_igblast_query(file_path=file_path)
    #     elif chain.lower() == 'heavy':
    #         self._heavy_chains.load_igblast_query(file_path=file_path)
    #     else:
    #         raise ValueError('Specify if the data being loaded is for the heavy or light chain')

    def numbering_table(self, as_array=False, region='all', chain='both'):

        return to_numbering_table(as_array=as_array, region=region, chain=chain,
                                  heavy_chains_numbering_table=self._heavy_chain.ab_numbering_table,
                                  light_chains_numbering_table=self._light_chain.ab_numbering_table,
                                  names=[self.name])

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._pair_sequence

    @property
    def germline_identity(self):
        return self._germline_identity()

    # @property
    # def germline(self):
    #
    #     heavy_germline = self._heavy_chains.germline
    #     light_germline = self._light_chains.germline
    #
    #     return {name_i: {"Heavy": heavy_germline[heavy_i],
    #                      "Light": light_germline[light_i]} for name_i, heavy_i,
    #                                                            light_i in zip(self._names, self._internal_heavy_name,
    #                                                                           self._internal_light_name)}

    def _string_summary_basic(self):
        return "abpytools.Fab Name: {} Sequence length: {}".format(self.name, len(self.sequence))

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __getitem__(self, index):
        if index == 0:
            # get light chain
            return self._light_chain
        elif index == 1:
            return self._heavy_chain
        else:
            raise ValueError("Can only slice object with either 0 (light chain) or 1 (heavy chain)")

    def _germline_identity(self):

        # empty dictionaries return false
        if bool(self[1].germline_identity) is False:
            # this means there is no information about the germline,
            # by default it will run a web query
            # this is a very lazy fix to to do a web query using a Chain object...
            ChainCollection(antibody_objects=[self[1]]).igblast_server_query()
        if bool(self[0].germline_identity) is False:
            ChainCollection(antibody_objects=[self[0]]).igblast_server_query()

        return germline_identity_pd({self._internal_heavy_name: self[0].germline_identity},
                                    {self._internal_light_name: self[1].germline_identity},
                                    [self._internal_heavy_name],
                                    [self._internal_light_name],
                                    [self._name])
