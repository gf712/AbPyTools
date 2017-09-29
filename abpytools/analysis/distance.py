from ..core.chain_collection import ChainCollection
import seaborn as sns
from .analysis_helper_functions import switch_interactive_mode
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import scipy.spatial.distance as ssd
from ..utils.python_config import PythonConfig


class DistancePlot(ChainCollection):
    def __init__(self, antibody_objects=None, path=None):
        super().__init__(antibody_objects=antibody_objects, path=path)

    def plot_heatmap(self, feature='chou', distance_metric='cosine_distance', save=False, ax=None, labels=None,
                     multiprocessing=False, **kwargs):

        data = self.distance_matrix(feature=feature, metric=distance_metric, multiprocessing=multiprocessing)

        switch_interactive_mode(save=save)

        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.set(xlabel='Antibody', ylabel='Antibody', title=distance_metric,
                   xticks=range(self.n_ab), yticks=range(self.n_ab))

        if labels is None:
            labels = self.names

        sns.heatmap(data, ax=ax, **kwargs)
        ax.set_yticklabels(labels, rotation='horizontal')
        ax.set_xticklabels(labels, rotation=60)

        ipython_config = PythonConfig()
        if ipython_config.ipython_info == 'notebook' and save is False:
            plt.plot()

    def plot_dendrogram(self, feature='chou', distance_metric='cosine_distance', save=False, ax=None, labels=None,
                        multiprocessing=False, **kwargs):

        switch_interactive_mode(save=save)

        data = self.distance_matrix(feature=feature, metric=distance_metric, multiprocessing=multiprocessing)
        # convert the redundant n*n square matrix form into a condensed nC2 array
        data = ssd.squareform(data)

        clustered_data = linkage(y=data)

        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.set(xlabel='Antibody', ylabel='Distance', title=distance_metric)

        if labels is None:
            labels = self.names

        # plot dendrogram
        _ = dendrogram(clustered_data, labels=labels, ax=ax, **kwargs)

        ipython_config = PythonConfig()
        if ipython_config.ipython_info == 'notebook' and save is False:
            plt.plot()
