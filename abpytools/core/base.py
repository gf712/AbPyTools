from abc import ABC, abstractmethod


class CollectionBase(ABC):

    """
    Abstract class to not forget some essential methods
    """

    @staticmethod
    @abstractmethod
    def load_from_file(path, n_threads, verbose, show_progressbar, **kwargs):
        pass  # pragma: no cover

    @abstractmethod
    def save(self, file_format, path, information):
        pass  # pragma: no cover
