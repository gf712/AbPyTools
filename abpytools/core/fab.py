from .chain_collection import Chain
import numpy as np
import pandas as pd
from .chain import calculate_charge
from abpytools.utils import DataLoader


class Fab:
    """
    
    """
    def __init__(self, heavy_chain=None, light_chain=None, load=True, name=None):
        # As a convention the order will always be light chain (index 0) and then heavy chain (index 1)

        # check if it's an Chain class
        if heavy_chain is None and light_chain is None:
            raise IOError('Provide a list of Chain objects or an ChainCollection object')

        if isinstance(heavy_chain, Chain):
            self._heavy_chain = heavy_chain

        else:
            raise ValueError('Provide a list of Chain objects or an ChainCollection object')

        if isinstance(light_chain, Chain):
            self._light_chain = light_chain

        else:
            raise ValueError('Provide a list of Chain objects or an ChainCollection object')

        # load data
        if load:
            self._heavy_chain.load()
            self._light_chain.load()

        # self._pair_sequence = self._heavy_chain.sequence + self._light_chain.sequence
        self._pair_sequence = self[0].sequence + self[1].sequence

        if isinstance(name, str):
            self._name = name

        elif name is None:
            self._name = 'ID1'

        # keep the name of the heavy and light chains internally to keep everything in the right order
        self._internal_heavy_name = self[1].name
        self._internal_light_name = self[0].name

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

    # def load_imgt_query(self, file_path, chain):
    #
    #     if chain.lower() == 'light':
    #         self._light_chains.load_imgt_query(file_path=file_path)
    #     elif chain.lower() == 'heavy':
    #         self._heavy_chains.load_imgt_query(file_path=file_path)
    #     else:
    #         raise ValueError('Specify if the data being loaded is for the heavy or light chain')

    def numbering_table(self, as_array=False, region='all', chain='both'):

        if chain == 'both':

            if as_array:
                t_heavy = self[1].ab_numbering_table(as_array=True, region=region)
                t_light = self[0].ab_numbering_table(as_array=True, region=region)

                data = np.concatenate((t_light, t_heavy), axis=1)

            else:
                t_heavy = self[1].ab_numbering_table(as_array=False, region=region)
                t_light = self[0].ab_numbering_table(as_array=False, region=region)

                t_heavy.reset_index(drop=True, inplace=True)
                t_light.reset_index(drop=True, inplace=True)

                data = pd.concat([t_light, t_heavy], axis=1, keys=['Light', 'Heavy'])

        elif chain == 'heavy':

            data = self[1].ab_numbering_table(as_array=as_array, region=region)

        elif chain == 'light':

            data = self[0].ab_numbering_table(as_array=as_array, region=region)

        else:
            raise ValueError("Unknown chain.")

        if not as_array:
            data.index = [self._name]

        return data

    # @property
    # def regions(self):
    #     heavy_regions = self._heavy_chain.ab_region_index()
    #     light_regions = self._light_chain.ab_region_index()
    #
    #     return {name: {'Heavy': heavy_regions[heavy], 'Light': light_regions[light]} for name, heavy, light in
    #             zip(self.names, self._internal_heavy_name, self._internal_light_name)}

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._pair_sequence

    # @property
    # def germline_identity(self):
    #     return self._germline_identity()
    #
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

        h_germline_pd = pd.DataFrame(self[1].germline_identity).T
        l_germline_pd = pd.DataFrame(self[0].germline_identity).T

        l_columns = pd.MultiIndex.from_tuples([('Light', x) for x in l_germline_pd.columns], names=['Chain', 'Region'])
        h_columns = pd.MultiIndex.from_tuples([('Heavy', x) for x in h_germline_pd.columns], names=['Chain', 'Region'])
        average_columns = pd.MultiIndex.from_tuples([('Average', x) for x in l_germline_pd.columns],
                                                    names=['Chain', 'Region'])

        l = pd.DataFrame(index=self._internal_light_name,
                         columns=l_columns)
        h = pd.DataFrame(index=self._internal_heavy_name,
                         columns=h_columns)

        l = l.apply(lambda x: l_germline_pd.ix[x.name], axis=1)
        h = h.apply(lambda x: h_germline_pd.ix[x.name], axis=1)

        average = (h.as_matrix() + l.as_matrix()) / 2

        l.columns = l_columns
        h.columns = h_columns
        average = pd.DataFrame(average, columns=average_columns, index=self._names)

        l.index = self._name
        h.index = self._name

        return pd.concat([h, l, average], axis=1)