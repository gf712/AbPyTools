from matplotlib import pyplot as plt
import seaborn as sns
from abpytools import AntibodyCollection
from os.path import join
from abpytools.utils import PythonConfig


class CDRLength:

    def __init__(self, path=None, antibody_collection=None, plot_path='./', plot_name='CDR_length',
                 plot_title='CDR Lengths', hist=True, notebook_plot=True, only_CDR3=True):

        if path is None and antibody_collection is None:
            raise IOError("Nothing to work with, please provide a path to a FASTA file or a AntibodyCollection object")

        if isinstance(path, str):
            self._path = path
        elif isinstance(antibody_collection, AntibodyCollection):
            self.antibody_collection = antibody_collection
            self._path = path
        else:
            raise IOError("Please provide a valid input format for path (str) or obj_list (AntibodyCollection)")

        if isinstance(plot_name, str) and isinstance(plot_path, str):
            self._plot_path = plot_path
            self._plot_name = plot_name
        else:
            raise IOError("Expected a string for file name and a string of the path to save plot")
        self._hist = hist
        self._plot_title = plot_title
        self._notebook_plot = notebook_plot
        self._only_CDR3 = only_CDR3

    def plot_cdr(self):

        if self._path is not None:
            self.antibody_collection = AntibodyCollection(path=self._path)
            self.antibody_collection.load()

        cdrs = self.antibody_collection.cdr_lengths()

        if self._only_CDR3:
            plt.title('CDR3 Length', size=18)
            sns.distplot(cdrs[:, 2].astype(int), hist=self._hist)
            plt.ylabel('Density', size=14)
            plt.xlabel('CDR Length', size=14)
        else:
            f, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 6))
            ax = ax.ravel()
            plt.suptitle('CDR Length', size=20)
            for i, cdr in enumerate(['CDR 1', 'CDR 2', 'CDR 3']):

                ax[i].set_title(cdr, size=16)
                sns.distplot(cdrs[:, i].astype(int), hist=self._hist, ax=ax[i])

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
            plt.savefig(join(self._plot_path, self._plot_name))
