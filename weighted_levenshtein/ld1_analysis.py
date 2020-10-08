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
    print('TCR pairs with LD=1:', df.shape[0])
    return df


def format_df_ld1():
    # Construct dataframe with TCR of LD=1 with corresponding epitopes
    file = "../data/CDR3_epitopes_LD1.csv"
    if Path(file).is_file():
        df = pd.read_csv(file, sep=';')
    else:
        tcr = []
        with open("../data/TRb_CD8_n42675_uppercase_CDR3_unique.csv", "r") as f:
            for line in f:
                tcr.append(line.strip())

        # Load epitope data
        df = pd.read_csv('../data/TRb_CD8_n42675_uppercase.csv', sep=";")
        df = df[['CDR3', 'Epitope_peptide']]
        # Epitopes -> columns
        one_hot = pd.get_dummies(df['Epitope_peptide'], prefix='', prefix_sep='')
        df = df.drop('Epitope_peptide', axis=1)
        df = df.join(one_hot)
        # one row = one unique TCR
        df = df.groupby(df['CDR3']).sum()
        df[df != 0] = 1 # to binary

        # Load LD=1 TCR data
        df_ld1 = pd.read_csv("../data/ld1_tcr_wld.csv", sep=";")
        # get unique TCR with LD=1
        tcr = set(df_ld1[['CDR3_pair1', 'CDR3_pair2']].values.T.ravel())
        print(len(tcr))

        # filter dataframe with only TCR LD=1
        print(df.head())
        df = df[df.index.isin(tcr)]
        df.to_csv(file, sep=';')
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


def compute_distance_for_each_LD1pair(plot=False):
    ld1_dist_file = "../data/ld1_tcr_wld_distances.csv"

    if Path(ld1_dist_file).is_file():
        df_tcr_pairs = pd.read_csv(ld1_dist_file, sep=";")

    else:
        # for each TCR pairs with LD=1, compute distance of epitope binding
        df_tcr_pairs = get_ld1()
        df_epi = format_df_ld1()
        dist_cols = ["Hamming", "Jaccard"]
        print("Computing distance")
        df_tcr_pairs[dist_cols] = df_tcr_pairs.apply(lambda row: compute_distance(row, df_epi), axis=1)
        df_tcr_pairs.to_csv(ld1_dist_file, sep=';', index=False)

    if plot:
        dft_plot = df_tcr_pairs.iloc[:, 2:].melt("wld", var_name='Metric', value_name='Distance')
        sns.lineplot(data=dft_plot, x="wld", y="Distance",
                         hue='Metric')
        plt.xlabel('WLD')
        plt.title("Distance of Epitope binding profile of each TCR pair (LD=1)")
        plt.savefig('../data/plot/epitope_distance_LD1.pdf')
        plt.savefig('../data/plot/epitope_distance_LD1.png')


compute_distance_for_each_LD1pair(plot=True)