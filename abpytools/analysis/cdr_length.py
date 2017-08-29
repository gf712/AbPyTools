from matplotlib import pyplot as plt
import seaborn as sns
from abpytools import ChainCollection
import os
from abpytools.utils import PythonConfig
from abpytools.features.cdr import CDR


class CDRLength:

    def __init__(self, path=None, antibody_collection=None, plot_path='./', plot_name='CDR_length',
                 plot_title='CDR Lengths', hist=True, notebook_plot=True):

        if path is None and antibody_collection is None:
            raise IOError("Nothing to work with, please provide a path to a FASTA file or a ChainCollection object")

        if isinstance(path, str):
            self._path = path
        elif isinstance(antibody_collection, ChainCollection):
            self.antibody_collection = antibody_collection
            self._path = path
        else:
            raise IOError("Please provide a valid input format for path (str) or obj_list (ChainCollection)")

        if isinstance(plot_name, str) and isinstance(plot_path, str):
            self._plot_path = plot_path
            self._plot_name = plot_name
        else:
            raise IOError("Expected a string for file name and a string of the path to save plot")
        self._hist = hist
        self._plot_title = plot_title
        self._notebook_plot = notebook_plot

    def plot_cdr(self, only_cdr3=True):

        if self._path is not None:
            self.antibody_collection = ChainCollection(path=self._path)
            self.antibody_collection.load()

        cdrs = CDR(antibodies=self.antibody_collection)
        cdr_lengths = cdrs.cdr_length()

        if only_cdr3:
            plt.title('CDR3 Length', size=18)
            sns.distplot(cdr_lengths[:, 2], hist=self._hist)
            plt.ylabel('Density', size=14)
            plt.xlabel('CDR Length', size=14)
        else:
            f, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 6))
            ax = ax.ravel()
            plt.suptitle('CDR Length', size=20)
            for i, cdr in enumerate(['CDR 1', 'CDR 2', 'CDR 3']):

                ax[i].set_title(cdr, size=16)
                sns.distplot(cdr_lengths[:, i], hist=self._hist, ax=ax[i])

                if i == 1:
                    ax[i].set_ylabel('Density', size=16)

                if i == 2:
                    ax[i].set_xlabel('CDR Length', size=16)

            plt.tight_layout()
            plt.subplots_adjust(top=0.9)

        ipython_config = PythonConfig()
        ipython_config.get_ipython_info()
        if ipython_config.backend == 'notebook' and self._notebook_plot:
            plt.plot()
        else:
            plt.savefig(os.path.join(self._plot_path, self._plot_name))
            plt.close()
