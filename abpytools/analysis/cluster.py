from sklearn import cluster, decomposition
import numpy as np
from abpytools import AntibodyCollection
from matplotlib import pyplot as plt


class Cluster:
    def __init__(self, antibodies, metric='hydrophobicity', clustering_method='kmeans',
                 decomposition='PCA'):

        if isinstance(antibodies, AntibodyCollection):
            self.antibodies = antibodies
        elif isinstance(antibodies, str):
            self.antibodies = AntibodyCollection(path=antibodies)
        if self.antibodies.n_ab == 0:
            self.antibodies.load()

        self.metric = metric
        self.clustering_method = clustering_method
        self.decomposition_method = decomposition
        self.cluster_assignment = np.zeros(self.antibodies.n_ab, dtype=int)
        self.cluster_assignment_dict = dict()
        self._data = None

    def _collect_data(self):
        if self.metric == 'hydrophobicity':
            return self.antibodies.hydrophobicity_matrix()

    def cluster(self, n_components=0.95, n_clusters=3):

        if self.decomposition_method == 'PCA':
            decomposition_obj = decomposition.PCA(n_components)

        self._data = decomposition_obj.fit_transform(self._collect_data())

        if self.clustering_method == 'kmeans':
            clustering_obj = cluster.KMeans(n_clusters=n_clusters)

        self.cluster_assignment = clustering_obj.fit_predict(self._data)

        for i, antibody_obj in enumerate(self.antibodies.antibody_objects):

            assignment = 'Cluster_{}'.format(self.cluster_assignment[i])

            if assignment not in self.cluster_assignment_dict:
                self.cluster_assignment_dict[assignment] = list()

            self.cluster_assignment_dict[assignment].append(antibody_obj)

    def plot_cluster(self):

        if len(self.cluster_assignment_dict) == 0:
            self.cluster()

        color = iter(plt.get_cmap('Vega20').colors)

        plt.figure(figsize=(8, 8))

        for assignment in np.unique(self.cluster_assignment):
            c = next(color)

            plt.scatter(self._data[self.cluster_assignment == assignment, 0],
                        self._data[self.cluster_assignment == assignment, 1],
                        c=c, label='Cluster {}'.format(assignment))

            plt.legend(loc='best', prop={"size": 14})
