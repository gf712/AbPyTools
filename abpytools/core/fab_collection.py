from .chain_collection import ChainCollection
import numpy as np
import pandas as pd
from .chain import calculate_charge
from abpytools.utils import DataLoader
from operator import itemgetter
from .fab import Fab


class FabCollection:
    """

    """

    def __init__(self, fab=None, heavy_chains=None, light_chains=None, load=True, names=None):

        # check if it's a Chain object
        if heavy_chains is None and light_chains is None and fab is None:
            raise IOError('Provide a list of Chain objects or an ChainCollection object')

        # check if fab object is a list and if all object are abpytools.Fab objects
        if isinstance(fab, list) and all(isinstance(fab_i, Fab) for fab_i in fab):
            self._fab = fab
            self._light_chains = [x[0] for x in self._fab]
            self._heavy_chains = [x[1] for x in self._fab]

        if fab is None and (heavy_chains is not None and light_chains is not None):

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

            if isinstance(names, list):
                if len(names) == self._heavy_chains.n_ab:
                    self._names = names
                else:
                    raise ValueError(
                        'Length of name list must be the same as length of heavy_chains/light chains lists')

            if names is None:
                self._names = ['{} - {}'.format(heavy, light) for heavy, light in zip(self._heavy_chains.names,
                                                                                      self._light_chains.names)]

        self._n_ab = self._light_chains.n_ab
        self._pair_sequences = [heavy + light for heavy, light in zip(self._heavy_chains.sequences,
                                                                      self._light_chains.sequences)]

        # keep the name of the heavy and light chains internally to keep everything in the right order
        self._internal_heavy_name = self._heavy_chains.names
        self._internal_light_name = self._light_chains.names

    # even though it makes more sense to draw all these values from the base Fab objects this is much slower
    # whenever self._n_ab > 1 it makes more sense to use the self._heavy_chain and self._light_chain containers
    # in all the methods
    # in essence the abpytools.Fab object is just a representative building block that could in future just
    # cache data and would then represent a speed up in the calculations
    def molecular_weight(self, monoisotopic=False):

        return [heavy + light for heavy, light in zip(self._heavy_chains.molecular_weights(monoisotopic=monoisotopic),
                                                      self._light_chains.molecular_weights(monoisotopic=monoisotopic))]

    def extinction_coefficient(self, extinction_coefficient_database='Standard', reduced=False, normalise=False,
                               **kwargs):

        heavy_ec = self._heavy_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        light_ec = self._light_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)

        if normalise:
            return [(heavy + light) / mw for heavy, light, mw in
                    zip(heavy_ec, light_ec, self.molecular_weight(**kwargs))]
        else:
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

    def imgt_local_query(self, file_path, chain):

        if chain.lower() == 'light':
            self._light_chains.imgt_local_query(file_path=file_path)
        elif chain.lower() == 'heavy':
            self._heavy_chains.imgt_local_query(file_path=file_path)
        else:
            raise ValueError('Specify if the data being loaded is for the heavy or light chain')

    def imgt_server_query(self, **kwargs):
        self._light_chains.imgt_server_query(**kwargs)
        self._heavy_chains.imgt_server_query(**kwargs)

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

    def _string_summary_basic(self):
        return "abpytools.FabCollection Number of sequences: {}".format(self._n_ab)

    def __len__(self):
        return self._n_ab

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __getitem__(self, indices):
        if isinstance(indices, int):
            return Fab(heavy_chain=self._heavy_chains[indices],
                       light_chain=self._light_chains[indices],
                       name=self.names[indices], load=False)
        else:
            return FabCollection(heavy_chains=list(itemgetter(*indices)(self._heavy_chains)),
                                 light_chains=list(itemgetter(*indices)(self._light_chains)),
                                 names=list(itemgetter(*indices)(self._names)), load=False)
        # if isinstance(indices, int):
        #     indices = [indices]
        #
        # return FabCollection(heavy_chains=list(itemgetter(*indices)(self._heavy_chains)),
        #                      light_chains=list(itemgetter(*indices)(self._light_chains)),
        #                      names=list(itemgetter(*indices)(self._names)), load=False)

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
