from abc import ABC, abstractmethod


class CollectionBase(ABC):

    @abstractmethod
    def molecular_weights(self, monoisotopic):
        pass

    @abstractmethod
    def extinction_coefficients(self, extinction_coefficient_database,
                                reduced, normalise):
        pass

    @abstractmethod
    def hydrophobicity_matrix(self):
        pass

    @abstractmethod
    def total_charge(self):
        pass

    @abstractmethod
    def charge(self):
        pass

    @abstractmethod
    def load(self, n_threads, verbose, show_progressbar, **kwargs):
        pass

    @abstractmethod
    def save(self):
        pass
