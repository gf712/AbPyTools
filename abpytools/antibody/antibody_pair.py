from .antibody_collection import AntibodyCollection
import numpy as np
import pandas as pd
from .antibody import calculate_charge
from abpytools.utils import DataLoader


class AntibodyPair:
    def __init__(self, heavy_chains=None, light_chains=None, load=True, names=None):

        # check if it's an Antibody class
        if heavy_chains is None and light_chains is None:
            raise IOError('Provide a list of Antibody objects or an AntibodyCollection object')

        if isinstance(heavy_chains, list):
            self._heavy_chains = AntibodyCollection(antibody_objects=heavy_chains)

        elif isinstance(heavy_chains, AntibodyCollection):
            self._heavy_chains = heavy_chains

        else:
            raise IOError('Provide a list of Antibody objects or an AntibodyCollection object')

        if isinstance(light_chains, list):
            self._light_chains = AntibodyCollection(antibody_objects=light_chains)

        elif isinstance(light_chains, AntibodyCollection):
            self._light_chains = light_chains

        else:
            raise IOError('Provide a list of Antibody objects or an AntibodyCollection object')

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

    @property
    def numbering_table(self):

        return pd.DataFrame(np.concatenate((self._heavy_chains.numbering_table().as_matrix(),
                                            self._light_chains.numbering_table().as_matrix()), axis=1),
                            columns=self._heavy_chains.numbering_table().columns.append(
                            self._light_chains.numbering_table().columns), index=self.names)

    @property
    def regions(self):
        return {'Heavy': self._heavy_chains.ab_region_index(), 'Light': self._light_chains.ab_region_index()}

    @property
    def names(self):
        return self._names

    @property
    def sequences(self):
        return self._pair_sequences

    @property
    def n_ab(self):
        return self._n_ab
