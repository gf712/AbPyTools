from .antibody import Antibody
from .antibody_collection import AntibodyCollection


class AntibodyPair:
    def __init__(self, heavy_chains=None, light_chains=None, load=True):

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

        self.pair_sequences = [heavy + light for heavy, light in zip(self._heavy_chains.sequences(),
                                                                     self._light_chains.sequences())]

    def pairs_mw(self, monoisotopic=False):

        return [heavy + light for heavy, light in zip(self._heavy_chains.molecular_weights(monoisotopic=monoisotopic),
                                                      self._light_chains.molecular_weights(monoisotopic=monoisotopic))]

    def pairs_ec(self, extinction_coefficient_database='Standard', reduced=False):

        heavy_ec = self._heavy_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        light_ec = self._light_chains.extinction_coefficients(
            extinction_coefficient_database=extinction_coefficient_database,
            reduced=reduced)
        return [heavy + light for heavy, light in zip(heavy_ec, light_ec)]

    def names(self):
        return ['{} - {}'.format(heavy, light) for heavy, light in zip(self._heavy_chains.names(),
                                                                       self._light_chains.names())]

    def sequences(self):
        return [heavy + light for heavy, light in zip(self._heavy_chains.sequences(), self._light_chains.sequences())]
