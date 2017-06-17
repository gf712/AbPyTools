from .antibody_collection import ChainCollection
import numpy as np
import pandas as pd
from .antibody import calculate_charge
from abpytools.utils import DataLoader


class Fab:
    def __init__(self, heavy_chains=None, light_chains=None, load=True, names=None):

        # check if it's an Chain class
        if heavy_chains is None and light_chains is None:
            raise IOError('Provide a list of Chain objects or an ChainCollection object')

        if isinstance(heavy_chains, list):
            self._heavy_chains = ChainCollection(antibody_objects=heavy_chains)

        elif isinstance(heavy_chains, ChainCollection):
            self._heavy_chains = heavy_chains

        else:
            raise IOError('Provide a list of Chain objects or an ChainCollection object')

        if isinstance(light_chains, list):
            self._light_chains = ChainCollection(antibody_objects=light_chains)

        elif isinstance(light_chains, ChainCollection):
            self._light_chains = light_chains

        else:
            raise IOError('Provide a list of Chain objects or an ChainCollection object')

        if self._light_chains.n_ab != self._heavy_chains.n_ab:
            raise ValueError('Number of heavy chains must be the same of light chains')

        # load data
        if load:
            self._heavy_chains.load()
            self._light_chains.load()

        self._pair_sequences = [heavy + light for heavy, light in zip(self._heavy_chains.sequences,
                                                                      self._light_chains.sequences)]

        if isinstance(names, list):
            if len(names) == self._heavy_chains.n_ab:
                self._names = names
            else:
                raise ValueError('Length of name list must be the same as length of heavy_chains/light chains lists')

        if names is None:
            self._names = ['{} - {}'.format(heavy, light) for heavy, light in zip(self._heavy_chains.names,
                                                                                  self._light_chains.names)]

        self._n_ab = self._light_chains.n_ab

        # keep the name of the heavy and light chains internally to keep everything in the right order
        self._internal_heavy_name = self._heavy_chains.names
        self._internal_light_name = self._light_chains.names

    def molecular_weight(self, monoisotopic=False):

        return [heavy + light for heavy, light in zip(self._heavy_chains.molecular_weights(monoisotopic=monoisotopic),
                                                      self._light_chains.molecular_weights(monoisotopic=monoisotopic))]

    def extinction_coefficient(self, extinction_coefficient_database='Standard', reduced=False):

        heavy_ec = self._heavy_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        light_ec = self._light_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        return [heavy + light for heavy, light in zip(heavy_ec, light_ec)]

    def hydrophobicity_matrix(self):

        return np.column_stack((self._heavy_chains.hydrophobicity_matrix(), self._light_chains.hydrophobicity_matrix()))

    def charge(self):

        return np.column_stack((self._heavy_chains.charge, self._light_chains.charge))

    def total_charge(self, ph=7.4, pka_database='Wikipedia'):

        available_pi_databases = ["EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell", "Wikipedia", "Lehninger",
                                  "Grimsley"]
        assert pka_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pka_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties', data=['pI', pka_database])
        pka_data = data_loader.get_data()

        return [calculate_charge(sequence=seq, ph=ph, pka_values=pka_data) for seq in self.sequences]

    def load_imgt_query(self, file_path, chain):

        if chain.lower() == 'light':
            self._light_chains.load_imgt_query(file_path=file_path)
        elif chain.lower() == 'heavy':
            self._heavy_chains.load_imgt_query(file_path=file_path)
        else:
            raise ValueError('Specify if the data being loaded is for the heavy or light chain')

    def numbering_table(self, as_array=False, region='all', chain='both'):

        if chain == 'both':

            if as_array:
                t_heavy = self._heavy_chains.numbering_table(as_array=True, region=region)
                t_light = self._light_chains.numbering_table(as_array=True, region=region)

                data = np.concatenate((t_light, t_heavy), axis=1)

            else:
                t_heavy = self._heavy_chains.numbering_table(as_array=False, region=region)
                t_light = self._light_chains.numbering_table(as_array=False, region=region)

                t_heavy.reset_index(drop=True, inplace=True)
                t_light.reset_index(drop=True, inplace=True)

                data = pd.concat([t_light, t_heavy], axis=1, keys=['Light', 'Heavy'])

        elif chain == 'heavy':

            data = self._heavy_chains.numbering_table(as_array=as_array, region=region)

        elif chain == 'light':

            data = self._light_chains.numbering_table(as_array=as_array, region=region)

        else:
            raise ValueError("Unknown chain.")

        if not as_array:
            data.index = self._names

        return data

    @property
    def regions(self):
        heavy_regions = self._heavy_chains.ab_region_index()
        light_regions = self._light_chains.ab_region_index()

        return {name: {'Heavy': heavy_regions[heavy], 'Light': light_regions[light]} for name, heavy, light in
                zip(self.names, self._internal_heavy_name, self._internal_light_name)}

    @property
    def names(self):
        return self._names

    @property
    def sequences(self):
        return self._pair_sequences

    @property
    def n_ab(self):
        return self._n_ab

    @property
    def germline_identity(self):
        return self._germline_identity()

    @property
    def germline(self):

        heavy_germline = self._heavy_chains.germline
        light_germline = self._light_chains.germline

        return {name_i: {"Heavy": heavy_germline[heavy_i],
                         "Light": light_germline[light_i]} for name_i, heavy_i,
                                                               light_i in zip(self._names, self._internal_heavy_name,
                                                                              self._internal_light_name)}

    def _germline_identity(self):

        h_germline_pd = pd.DataFrame(self._heavy_chains.germline_identity).T
        l_germline_pd = pd.DataFrame(self._light_chains.germline_identity).T

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

        l.index = self._names
        h.index = self._names

        return pd.concat([h, l, average], axis=1)
