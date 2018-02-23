from abc import ABC, abstractmethod


class CollectionBase(ABC):

    """
    Abstract class to not forget some essential methods
    """

    @abstractmethod
    def load(self, n_threads, verbose, show_progressbar, **kwargs):
        pass  # pragma: no cover

    @abstractmethod
    def save(self, file_format, file_path, file_name, information):
        pass  # pragma: no cover
