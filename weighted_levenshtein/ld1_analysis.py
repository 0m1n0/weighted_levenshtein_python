from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import distance
from matplotlib import pyplot as plt
import seaborn as sns
pd.options.mode.chained_assignment = None  # default='warn'


def get_ld1(plot=False):
    # Returns TCR sequences with LD=1
    ld1_file = "../data/ld1_tcr_wld.csv"
    dist_ld_file = '../data/dist_ld.npy'
    dist_wld_file = '../data/dist_wld.npy'
    tcr_file = "../data/TRb_CD8_n42675_uppercase_CDR3_unique.csv"

    if Path(ld1_file).is_file():
        df = pd.read_csv(ld1_file, sep=";")

    else:
        print('Loading raw LD distance data')
        dist = np.load(dist_ld_file)
        print('Raw distance matrix shape:', dist.shape)
        # get pairwise index of matrix where value = 1
        idx = np.where(np.triu(dist) == 1)

        # get WLD of same index
        dist = np.load(dist_wld_file)
        wld = dist[idx]

        # load TCR data
        all_tcr = []
        with open(tcr_file, "r") as f:
            for line in f:
                all_tcr.append(line.strip())
        tcr1 = [all_tcr[i] for i in idx[0]]
        tcr2 = [all_tcr[i] for i in idx[1]]

        df = pd.DataFrame({'CDR3_pair1': tcr1, 'CDR3_pair2': tcr2, 'wld': wld})
        df.to_csv(ld1_file, sep=';', index=False)

    if plot:
        plt.hist(df['wld'], bins=40)
        plt.title("WLD distribution for LD=1")
        plt.xlabel("WLD")
        plt.ylabel("TCR pairs")
        plt.savefig('../data/plot/LD1_WLD_distribution.pdf')
        plt.savefig('../data/plot/LD1_WLD_distribution.png')
        plt.close()
    print('TCR pairs with LD=1:', df.shape[0])
    return df


def format_df_ld1(col):
    # Construct dataframe with TCR of LD=1 with corresponding epitopes
    file = "../data/CDR3_{}.csv".format(col)
    file_ld1 = "../data/CDR3_{}_LD1.csv".format(col)
    file_oneHot = "../data/CDR3_{}_LD1_oneHot.csv".format(col)
    if Path(file).is_file():
        print(file)
        df = pd.read_csv(file, sep=';')
    else:
        if col == 'Epitope':
            col = 'Epitope_peptide'
        tcr = []
        with open("../data/TRb_CD8_n42675_uppercase_CDR3_unique.csv", "r") as f:
            for line in f:
                tcr.append(line.strip())

        # Load epitope or pathology data
        df = pd.read_csv('../data/TRb_CD8_n42675_uppercase.csv', sep=";")
        df = df[['CDR3', col]]
        # Epitopes or Pathology -> columns
        one_hot = pd.get_dummies(df[col], prefix='', prefix_sep='')
        df = df.drop(col, axis=1)
        df = df.join(one_hot)
        # one row = one unique TCR
        df = df.groupby(df['CDR3']).sum()
        df.to_csv(file, sep=';')

        # Load LD=1 TCR data
        df_ld1 = pd.read_csv("../data/ld1_tcr_wld.csv", sep=";")
        # get unique TCR with LD=1
        tcr = set(df_ld1[['CDR3_pair1', 'CDR3_pair2']].values.T.ravel())
        print(len(tcr))

        # filter dataframe with only TCR LD=1
        print(df.head())
        df = df[df.index.isin(tcr)]
        df.to_csv(file_ld1, sep=';')
        df[df != 0] = 1  # to binary
        df.to_csv(file_oneHot, sep=';')
    print("Unique TCR:", df.shape[0])
    return df


def compute_distance(row, df_epi):
    # return Hamming and Jaccard distance
    tcr_pair = [row['CDR3_pair1'], row['CDR3_pair2']]
    df_epi_pair = df_epi[df_epi['CDR3'].isin(tcr_pair)].iloc[:, 1:]
    vec1, vec2 = df_epi_pair.values.tolist()[0], df_epi_pair.values.tolist()[1]
    hamming = distance.hamming(vec1, vec2)
    jaccard = distance.jaccard(vec1, vec2)
    return pd.Series((hamming, jaccard))


def compute_distance_for_each_LD1pair(col, plot=False):
    ld1_dist_file = "../data/ld1_tcr_wld_distances_{}.csv".format(col)

    if Path(ld1_dist_file).is_file():
        df_tcr_pairs = pd.read_csv(ld1_dist_file, sep=";")

    else:
        # for each TCR pairs with LD=1, compute distance of epitope binding or pathology
        df_tcr_pairs = get_ld1()
        df_epi = format_df_ld1(col)
        dist_cols = ["Hamming", "Jaccard"]
        print("Computing distance")
        df_tcr_pairs[dist_cols] = df_tcr_pairs.apply(lambda row: compute_distance(row, df_epi), axis=1)
        df_tcr_pairs.to_csv(ld1_dist_file, sep=';', index=False)

    if plot:
        dft_plot = df_tcr_pairs.iloc[:, 2:].melt("wld", var_name='Metric', value_name='Distance')
        sns.lineplot(data=dft_plot, x="wld", y="Distance",
                         hue='Metric')
        plt.xlabel('WLD')
        plt.title("Distance of {} binding profile of each TCR pair (LD=1)".format(col))
        plt.savefig('../data/plot/distance_LD1_{}.pdf'.format(col))
        plt.savefig('../data/plot/distance_LD1_{}.png'.format(col))
        plt.close()


get_ld1(plot=True)
# col='Epitope' or 'Pathology'
compute_distance_for_each_LD1pair(col='Epitope', plot=True)
compute_distance_for_each_LD1pair(col='Pathology', plot=True)