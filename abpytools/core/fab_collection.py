from .chain_collection import ChainCollection
import numpy as np
import pandas as pd
from .chain import calculate_charge
from abpytools.utils import DataLoader
from operator import itemgetter
from .fab import Fab
from .helper_functions import germline_identity_pd, to_numbering_table
from .base import CollectionBase


class FabCollection(CollectionBase):
    """

    """

    def __init__(self, fab=None, heavy_chains=None, light_chains=None, names=None):

        # check if it's a Chain object
        if heavy_chains is None and light_chains is None and fab is None:
            raise ValueError('Provide a list of Chain objects or an ChainCollection object')

        # check if fab object is a list and if all object are abpytools.Fab objects
        if isinstance(fab, list) and all(isinstance(fab_i, Fab) for fab_i in fab):
            self._fab = fab
            self._light_chains = ChainCollection([x[0] for x in self._fab])
            self._heavy_chains = ChainCollection([x[1] for x in self._fab])

        if fab is None and (heavy_chains is not None and light_chains is not None):

            if isinstance(heavy_chains, list):
                self._heavy_chains = ChainCollection(antibody_objects=heavy_chains)

            elif isinstance(heavy_chains, ChainCollection):
                self._heavy_chains = heavy_chains

            else:
                raise ValueError('Provide a list of Chain objects or an ChainCollection object')

            if isinstance(light_chains, list):
                self._light_chains = ChainCollection(antibody_objects=light_chains)

            elif isinstance(light_chains, ChainCollection):
                self._light_chains = light_chains

            else:
                raise ValueError('Provide a list of Chain objects or an ChainCollection object')

            if len(self._light_chains.loading_status()) == 0:
                self._light_chains.load()

            if len(self._heavy_chains.loading_status()) == 0:
                self._heavy_chains.load()

            if self._light_chains.n_ab != self._heavy_chains.n_ab:
                raise ValueError('Number of heavy chains must be the same of light chains')

        if isinstance(names, list) and all(isinstance(name, str) for name in names):
            if len(names) == self._heavy_chains.n_ab:
                self._names = names
            else:
                raise ValueError(
                    'Length of name list must be the same as length of heavy_chains/light chains lists')

        elif names is None:
            self._names = ['{} - {}'.format(heavy, light) for heavy, light in zip(self._heavy_chains.names,
                                                                                  self._light_chains.names)]

        else:
            raise ValueError("Names expected a list of strings, instead got {}".format(type(names)))

        self._n_ab = self._light_chains.n_ab
        self._pair_sequences = [heavy + light for light, heavy in zip(self._heavy_chains.sequences,
                                                                      self._light_chains.sequences)]

        # keep the name of the heavy and light chains internally to keep everything in the right order
        self._internal_heavy_name = self._heavy_chains.names
        self._internal_light_name = self._light_chains.names

    # even though it makes more sense to draw all these values from the base Fab objects this is much slower
    # whenever self._n_ab > 1 it makes more sense to use the self._heavy_chain and self._light_chain containers
    # in all the methods
    # in essence the abpytools.Fab object is just a representative building block that could in future just
    # cache data and would then represent a speed up in the calculations
    def molecular_weights(self, monoisotopic=False):

        return [heavy + light for heavy, light in zip(self._heavy_chains.molecular_weights(monoisotopic=monoisotopic),
                                                      self._light_chains.molecular_weights(monoisotopic=monoisotopic))]

    def extinction_coefficients(self, extinction_coefficient_database='Standard', reduced=False, normalise=False,
                                **kwargs):

        heavy_ec = self._heavy_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        light_ec = self._light_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)

        if normalise:
            return [(heavy + light) / mw for heavy, light, mw in
                    zip(heavy_ec, light_ec, self.molecular_weights(**kwargs))]
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

    def igblast_local_query(self, file_path, chain):

        if chain.lower() == 'light':
            self._light_chains.igblast_local_query(file_path=file_path)
        elif chain.lower() == 'heavy':
            self._heavy_chains.igblast_local_query(file_path=file_path)
        else:
            raise ValueError('Specify if the data being loaded is for the heavy or light chain')

    def igblast_server_query(self, **kwargs):
        self._light_chains.igblast_server_query(**kwargs)
        self._heavy_chains.igblast_server_query(**kwargs)

    def numbering_table(self, as_array=False, region='all', chain='both', **kwargs):

        return to_numbering_table(as_array=as_array, region=region, chain=chain,
                                  heavy_chains_numbering_table=self._heavy_chains.numbering_table,
                                  light_chains_numbering_table=self._light_chains.numbering_table,
                                  names=self.names, **kwargs)

    def _germline_pd(self):

        heavy_chain_germlines = self._heavy_chains.germline
        light_chain_germlines = self._light_chains.germline

        data = np.array([[heavy_chain_germlines[x][0] for x in self._internal_heavy_name],
                         [heavy_chain_germlines[x][1] for x in self._internal_heavy_name],
                         [light_chain_germlines[x][0] for x in self._internal_light_name],
                         [light_chain_germlines[x][1] for x in self._internal_light_name]]).T

        df = pd.DataFrame(data=data,
                          columns=pd.MultiIndex.from_tuples([('Heavy', 'Assignment'),
                                                             ('Heavy', 'Score'),
                                                             ('Light', 'Assignment'),
                                                             ('Light', 'Score')]),
                          index=self.names)

        df = df.loc[:, (slice(None), 'Score')].apply(pd.to_numeric)

        return df

    def save(self):
        # TODO
        raise NotImplementedError("This is not the code you are looking for.")

    def load(self):
        # TODO
        raise NotImplementedError("This is not the code you are looking for.")

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
        return self._germline_pd()

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
                                 names=list(itemgetter(*indices)(self._names)))
        # if isinstance(indices, int):
        #     indices = [indices]
        #
        # return FabCollection(heavy_chains=list(itemgetter(*indices)(self._heavy_chains)),
        #                      light_chains=list(itemgetter(*indices)(self._light_chains)),
        #                      names=list(itemgetter(*indices)(self._names)), load=False)

    def _germline_identity(self):

        # empty dictionaries return false, so this condition checks if any of the values are False
        if all([x for x in self._light_chains.germline_identity.values()]) is False:
            # this means there is no information about the germline,
            # by default it will run a web query
            self._light_chains.igblast_server_query()
        if all([x for x in self._heavy_chains.germline_identity.values()]) is False:
            self._heavy_chains.igblast_server_query()

        return germline_identity_pd(self._heavy_chains.germline_identity,
                                    self._light_chains.germline_identity,
                                    self._internal_heavy_name,
                                    self._internal_light_name,
                                    self._names)

    def get_object(self, name):

        """

        :param name: str
        :return:
        """

        if name in self.names:
            index = self.names.index(name)
            return self[index]
        else:
            raise ValueError('Could not find sequence with specified name')
