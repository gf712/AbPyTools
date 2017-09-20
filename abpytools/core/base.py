from abc import ABC, abstractmethod


class CollectionBase(ABC):

    @abstractmethod
    def molecular_weights(self):
        pass

    @abstractmethod
    def extinction_coefficients(self):
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
    def load(self):
        pass

    @abstractmethod
    def save(self):
        pass
