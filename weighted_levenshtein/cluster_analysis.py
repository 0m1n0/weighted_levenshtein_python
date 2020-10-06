import pandas as pd
from pathlib import Path


def tcr_cluster_df(dist_type):
    filtered_tcr_file = "../data/filtered_tcr_{}.txt".format(dist_type)
    tcr = []
    with open(filtered_tcr_file, "r") as f:
        for line in f:
            tcr.append(line.strip())

    cluster_file = '../data/dbscan_clusters_{}.txt'.format(dist_type)
    clusters = []
    with open(cluster_file, "r") as f:
        for line in f:
            clusters.append(int(line.strip()))

    df = pd.DataFrame({'CDR3': tcr,
                       'cluster': clusters})
    df = df.set_index('CDR3')
    return df


def merge_df():
    file = "../data/CDR3_cluster_epitopes.csv"
    if Path(file).is_file():
        df = pd.read_csv(file, sep=';')
    else:
        # join LD & WLD clusters data
        dfl = tcr_cluster_df('ld')
        dfw = tcr_cluster_df('wld')
        df_cluster = dfl.join(dfw, lsuffix='_LD', rsuffix='_WLD')

        # Load epitope data
        df = pd.read_csv('../data/TRb_CD8_n42675_uppercase.csv', sep=";")
        df = df[['CDR3', 'Epitope_peptide']]
        # Epitopes -> columns
        one_hot = pd.get_dummies(df['Epitope_peptide'], prefix='', prefix_sep='')
        df = df.drop('Epitope_peptide', axis=1)
        df = df.join(one_hot)
        # CDR3 -> index
        df = df.set_index('CDR3')
        # keep unique CDR3
        df = df.groupby(df.index).sum()
        df[df != 0] = 1
        df['epi_sum'] = df.sum(axis=1)

        df = df_cluster.join(df)
        df = df.sort_values(by=['cluster_LD'])

        df.to_csv(file, sep=';')
    return df

merge_df()

# https://scikit-learn.org/stable/modules/model_evaluation.html#multilabel-ranking-metrics
# 3.3.3.1. Coverage error
import numpy as np
from sklearn.metrics import coverage_error
y_true = np.array([[1, 0, 0], [0, 0, 1]])

y_score = np.array([[0.75, 0.5, 1], [1, 0.2, 0.1]])
coverage_error(y_true, y_score)
2.5

y_score = np.array([[1, 0, 0], [0, 0, 1]])
1

y_score = np.array([[4, 0, 0], [0, 0, 6]])
1

y_score = np.array([[0.4, 0, 0], [0, 0, 60]])
1

y_score = np.array([[0, 0, 0], [1, 0, 0]])
3

y_score = np.array([[0, 0, 0], [100, 0, 0]])
3

y_true = np.array([[1, 0, 1], [0, 0, 1]])
y_score = np.array([[4, 0, 1], [0, 0, 1]])
1.5

y_score = np.array([[4, 1, 1], [0, 0, 1]])
2

# a[a['CDR3'].duplicated()]
#              CDR3  AARAVFLAL  GILGFVFTL  ...  RPRGEVRFL  RYPLTFGW  RYPLTFGWCF
# 9   CASSIRSSYEQYF          0          1  ...          0         0           0
# 11  CASSIRSSYEQYF          0          1  ...          0         0           0
# 56  CASSSRSSYEQYF          0          1  ...          0         0           0
# 64  CASSIRSTDTQYF          0          1  ...          0         0           0
#
#
# CSARGQRGLVNEQFF
#
# b = a[:12]
# b = b.set_index('CDR3')
# b = b.groupby(b.index).sum()
# b[b != 0] = 1
# b['sum'] = b.sum(axis=1)
#
