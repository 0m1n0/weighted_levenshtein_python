from pathlib import Path
import numpy as np
import pandas as pd
import csv
from sklearn.datasets.samples_generator import make_blobs
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()


def distance_matrix(dist_type):
    filtered_dist_file = "../data/filtered_dist_{}.npy".format(dist_type)
    filtered_tcr_file = "../data/filtered_tcr_{}.txt".format(dist_type)

    if Path(filtered_dist_file).is_file():
        print('Loading filtered distance data')
        dist = np.load(filtered_dist_file)
        tcr = []
        with open(filtered_tcr_file, "r") as f:
            for line in f:
                tcr.append(line.strip())
        print('Filtered distance matrix shape:', dist.shape)

    else:
        # Load Weighted or not Levenshtein distance matrix
        print('Loading raw distance data')
        dist = np.load('../data/dist_{}.npy'.format(dist_type))
        print('Raw distance matrix shape:', dist.shape)

        # Load TCR sequence data
        with open('../data/TRb_CD8_n42675_uppercase_CDR3.csv', 'r') as f:
            tcr = [j for i in list(csv.reader(f)) for j in i]

        # max weighted Levenshtein distance of a substitution is 15
        min_value = 1 if dist_type == "ld" else 15
        print("Distance threshold:", min_value)
        # Select rows with minimum value
        rows = np.where(np.any((dist > 0) & (dist <= min_value), axis=1))[0]
        print("Removed fraction:", 1 - rows.shape[0]/dist.shape[0])
        dist = dist[np.ix_(rows, rows)]  # symmetric
        print('Filtered distance matrix shape:', dist.shape)


        # Get TCR corresponding to selected rows
        tcr = [tcr[i] for i in rows]

        # Save
        print('Saving')
        np.save(filtered_dist_file, dist)
        with open(filtered_tcr_file, "w") as f:
            for i in tcr:
                f.write(str(i) + "\n")
    return dist, tcr


def remove_duplicates(dist_type, dist, tcr):
    filtered_dist_file = "../data/filtered_dist_{}.npy".format(dist_type)
    filtered_tcr_file = "../data/filtered_tcr_{}.txt".format(dist_type)

    print('Distance matrix shape:', dist.shape)
    unique_idx = pd.Series(tcr).drop_duplicates().index
    tcr = pd.Series(tcr).drop_duplicates().to_list()
    dist = dist[np.ix_(unique_idx, unique_idx)]
    print('After drop of duplicated TCR sequences:', dist.shape)

    # Save
    print('Saving')
    np.save(filtered_dist_file, dist)
    with open(filtered_tcr_file, "w") as f:
        for i in tcr:
            f.write(str(i) + "\n")

    return dist, tcr


def get_optimal_epsilon(dist_type, dist):
    plot_file = '../data/plot/NearestNeighbors_{}.png'.format(dist_type)
    if Path(plot_file).is_file():
        print("Nearest Neighbors already computed")

    else:
        # Compute optimal value for epsilon of DBSCAN
        # Nearest Neighbors algorighm used
        # https://towardsdatascience.com/machine-learning-clustering-dbscan-determine-the-optimal-value-for-epsilon-eps-python-example-3100091cfbc
        print("Using KNN to find optimal value for epsilon of DBSCAN")
        neigh = NearestNeighbors(n_neighbors=2, metric='euclidean')
        nbrs = neigh.fit(dist)
        # distance to the closest n-neighbors points
        distances, indices = nbrs.kneighbors(dist)
        # sort and plot results
        distances = np.sort(distances, axis=0)
        distances = distances[:, 1]
        plt.plot(distances)
        plt.savefig(plot_file)
        plt.savefig('../data/plot/NearestNeighbors_{}.pdf'.format(dist_type))



def dbscan(dist_type, dist):
    # epsilon = 2 if dist_type == "ld" else 16
    # optimal epsilon between 0.1-40% of lowest distance -> mean
    sort = np.unique(dist)
    epsilon = round((sort[int(len(sort) * 0.01)] + sort[int(len(sort) * 0.4)]) / 2)
    print('epsilon value =', epsilon)
    clust = DBSCAN(eps=epsilon, metric="precomputed", min_samples=2)
    clust.fit(dist)
    # list of clusters and their respective points
    clusters = clust.labels_
    np.savetxt('../data/dbscan_clusters_{}.txt'.format(dist_type), clusters,
               delimiter=',', fmt='%i')


dist_type = "ld"
dist, tcr = distance_matrix(dist_type)
dist, tcr = remove_duplicates(dist_type, dist, tcr)
dbscan(dist_type, dist)

dist_type = "wld"
dist, tcr = distance_matrix(dist_type)
dist, tcr = remove_duplicates(dist_type, dist, tcr)
dbscan(dist_type, dist)