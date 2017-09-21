from matplotlib import pyplot as plt
import seaborn as sns
import os
from abpytools.utils import PythonConfig
from abpytools.features.regions import ChainDomains
from matplotlib.ticker import MaxNLocator


class CDRLength(ChainDomains):

    def __init__(self, path=None, antibody_objects=None, load=True):

        super().__init__(path=path, antibody_objects=antibody_objects, load=load)

    def plot_cdr(self, only_cdr3=True, save=False, plot_path='./', plot_name='CDR_length',
                 plot_title=None, hist=True, ax=None, **kwargs):

        ipython_config = PythonConfig()
        if ipython_config.ipython_info == 'notebook' and save is False:
            if ipython_config.matplotlib_interactive is False:
                # turns on interactive mode
                plt.ion()

        if ax is None:
            if only_cdr3:
                f, ax = plt.subplots(nrows=1, ncols=1)
            else:
                f, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
                ax = ax.ravel()

        if only_cdr3:
            if plot_title is None:
                ax.set_title('CDR3 Length', size=18)
            else:
                ax.set_title(plot_title, size=18)
            sns.distplot(self.cdr_lengths()[:, 2], hist=hist, ax=ax, **kwargs)
            ax.set_ylabel('Density', size=14)
            ax.set_xlabel('CDR Length', size=14)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            if plot_title is None:
                plt.suptitle('CDR Length', size=20)
            else:
                plt.suptitle(plot_title, size=20)
            for i, cdr in enumerate(['CDR 1', 'CDR 2', 'CDR 3']):

                ax[i].set_title(cdr, size=16)
                sns.distplot(self.cdr_lengths()[:, i], hist=hist, ax=ax[i])

                if i == 0:
                    ax[i].set_ylabel('Density', size=16)

                if i == 1:
                    ax[i].set_xlabel('CDR Length', size=16)

                ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

            plt.tight_layout()
            plt.subplots_adjust(top=0.85)

        if ipython_config.ipython_info == 'notebook' and save is False:
            plt.plot()
        else:
            plt.savefig(os.path.join(plot_path, plot_name), format='png')
            plt.close()
