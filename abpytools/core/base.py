from abc import ABC, abstractmethod
import os
from .flags import *


class CollectionBase(ABC):

    """
    CollectionBase is the abpytools base class to develop its API
    """
    @classmethod
    def load_from_json(cls, path, n_threads=20, verbose=True, show_progressbar=True):
        raise NotImplementedError

    @classmethod
    def load_from_pb2(cls, path, n_threads=20, verbose=True, show_progressbar=True):
        raise NotImplementedError

    @classmethod
    def load_from_fasta(cls, path, numbering_scheme=NUMBERING_FLAGS.CHOTHIA, n_threads=20,
                        verbose=True, show_progressbar=True):
        raise NotImplementedError

    @classmethod
    def load_from_file(cls, path, n_threads=20, verbose=True, show_progressbar=True, **kwargs):

        """
        Args:
            path:
            n_threads: int to specify number of threads to use in loading process
            verbose: bool controls the level of verbose
            show_progressbar: bool whether to display the progressbar
            kwargs:

        Returns:
              ChainCollection
        """
        # check if path to file is valid
        if not os.path.isfile(path):
            raise ValueError("File does not exist!")

        file_format = path.split('.')[-1]
        if file_format not in ['json', 'pb2', 'fasta']:
            raise ValueError("Expected the file format to be json, pb2 or fasta.")

        if file_format == "json":
            chain_collection = cls.load_from_json(path,
                                                  n_threads=n_threads,
                                                  verbose=verbose,
                                                  show_progressbar=show_progressbar)

        elif file_format == "pb2" and BACKEND_FLAGS.HAS_PROTO:
            chain_collection = cls.load_from_pb2(path,
                                                 n_threads=n_threads,
                                                 verbose=verbose,
                                                 show_progressbar=show_progressbar)

        elif file_format == "pb2" and not BACKEND_FLAGS.HAS_PROTO:
            raise ValueError("protobuf has to be enabled to serialise objects with protobuf")

        elif file_format == "fasta":
            chain_collection = cls.load_from_fasta(path,
                                                   n_threads=n_threads,
                                                   verbose=verbose,
                                                   show_progressbar=show_progressbar,
                                                   **kwargs)
        else:
            raise ValueError("Expected the file format to be json, pb2 or fasta.")

        return chain_collection

    @abstractmethod
    def save(self, file_format, path, information):
        pass  # pragma: no cover
