from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


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

    print('TCR with LD=1:', df.shape[0])
    return df


df = get_ld1(plot=True)
print(df['wld'].describe())

