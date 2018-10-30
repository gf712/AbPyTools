# from abc import ABC, abstractmethod
import os
from .flags import *

AVAILABLE_FORMATS = [FORMAT_FLAGS.JSON, FORMAT_FLAGS.FASTA, FORMAT_FLAGS.PB2]


class CollectionBase:
    """
    CollectionBase is the abpytools base class to develop the collection APIs
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

        """
        # check if path to file is valid
        if not os.path.isfile(path):
            raise ValueError("File does not exist!")

        file_format = path.split('.')[-1]
        if file_format not in AVAILABLE_FORMATS:
            raise ValueError("Expected the file format to be json, pb2 or fasta.")

        if file_format == FORMAT_FLAGS.JSON:
            collection = cls.load_from_json(path,
                                            n_threads=n_threads,
                                            verbose=verbose,
                                            show_progressbar=show_progressbar)

        elif file_format == FORMAT_FLAGS.PB2 and BACKEND_FLAGS.HAS_PROTO:
            collection = cls.load_from_pb2(path,
                                           n_threads=n_threads,
                                           verbose=verbose,
                                           show_progressbar=show_progressbar)

        elif file_format == FORMAT_FLAGS.PB2 and not BACKEND_FLAGS.HAS_PROTO:
            raise ValueError("protobuf has to be enabled to serialise objects with protobuf")

        elif file_format == FORMAT_FLAGS.FASTA:
            collection = cls.load_from_fasta(path,
                                             n_threads=n_threads,
                                             verbose=verbose,
                                             show_progressbar=show_progressbar,
                                             **kwargs)
        else:
            raise NotImplementedError

        return collection

    def save_to_json(self, path, update=True):
        raise NotImplementedError

    def save_to_pb2(self, path, update=True):
        raise NotImplementedError

    def save_to_fasta(self, path, update=True):
        raise NotImplementedError

    def save(self, file_format, path, update=True):
        """

        Args:
            file_format:
            path:
            update:

        Returns:

        """

        # check if path is for a new or existing file
        if os.path.isfile(path) and update:
            # read old data and append to it
            update = True
        else:
            # overwrite to path
            update = False

        if file_format not in AVAILABLE_FORMATS:
            raise ValueError("Expected the file format to be json, pb2 or fasta.")

        if file_format == FORMAT_FLAGS.JSON:
            collection = self.save_to_json(path,
                                           update=update)

        elif file_format == FORMAT_FLAGS.PB2 and BACKEND_FLAGS.HAS_PROTO:
            collection = self.save_to_pb2(path,
                                          update=update)

        elif file_format == FORMAT_FLAGS.PB2 and not BACKEND_FLAGS.HAS_PROTO:
            raise ValueError("protobuf has to be enabled to serialise objects with protobuf")

        elif file_format == FORMAT_FLAGS.FASTA:
            collection = self.save_to_fasta(path,
                                            update=update)
        else:
            raise NotImplementedError

        return collection
