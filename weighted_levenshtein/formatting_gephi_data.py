import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
pd.options.mode.chained_assignment = None  # default='warn'


def correlation_plot(df, title):
    sns.set_theme(style="white")
    # Compute the correlation matrix
    corr = df.corr(method='spearman')
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))
    # Set up the matplotlib figure
    # f, ax = plt.subplots(figsize=(11, 9))
    plt.subplots(figsize=(11, 9))
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig('../data/plot/correlation_spearman_{}.pdf'.format(title))
    plt.savefig('../data/plot/correlation_spearman_{}.png'.format(title))


# formatting data for Gephi network visualization software
def rescale_wld_to_similarity(x):
    # for LD=1, we have WLD between 8-15
    wld_list = list(range(8, 16))
    # rescaled: between 1-8 as similarity score
    rescaled_value = x - 2 * wld_list.index(x)
    return rescaled_value


def create_edge_data(plot=False):
    file = "../data/gephi_LD1_edges.csv"
    if Path(file).is_file():
        print("Edge data already exists")
        df = pd.read_csv(file, sep=";")
    else:
        # TCR-TCR edges with weight
        df_tcr = pd.read_csv("../data/ld1_tcr_wld.csv", sep=";")
        # get unique TCR list
        tcr_list = (set(df_tcr['CDR3_pair1']).union(set(df_tcr['CDR3_pair2'])))
        df_tcr['rescaled_wld'] = df_tcr.apply(lambda row: rescale_wld_to_similarity(row['wld']), axis=1)
        df_tcr = pd.DataFrame({'Source': df_tcr['CDR3_pair1'],
                               'Target': df_tcr['CDR3_pair2'],
                               'Weight': df_tcr['rescaled_wld'],
                               'Edge_type': 'TCR-TCR'})
        print("Unique TCR edges:", df_tcr.shape[0])

        # TCR-Epitope edges with weight corresponding to pathology
        df_epi = pd.read_csv("../data/TRb_CD8_n42675_uppercase.csv", sep=";")
        # keep only TCR with LD=1
        df_epi = df_epi[df_epi['CDR3'].isin(tcr_list)]
        # get TCR-Pathology profile data
        df_tcr_patho = pd.read_csv("../data/CDR3_Pathology.csv", sep=";", index_col='CDR3')
        # add weight
        df_epi['Weight'] = df_epi.apply(lambda row: df_tcr_patho.loc[row['CDR3'], row['Pathology']], axis=1)

        df_epi = pd.DataFrame({'Source': df_epi['CDR3'],
                               'Target': df_epi['Epitope_peptide'],
                               'Weight': df_epi['Weight'],
                               'Edge_type': 'TCR-Epitope' })
        # concatenate TCR-Epitope
        df = pd.concat([df_tcr, df_epi], ignore_index=True)
        print("TCR + Epitope edges:", df.shape[0])
        df.to_csv(file, sep=';', index=False)

    if plot:
        ax = sns.countplot(x="Weight", hue="Edge_type", data=df)
        ax.set(yscale="log")
        plt.savefig('../data/plot/edge_weight_distribution.pdf')
        plt.savefig('../data/plot/edge_weight_distribution.png')
        plt.close()
    return df


def create_node_data(df_edge):
    file = "../data/gephi_LD1_nodes.csv"
    if Path(file).is_file():
        print("Node data already exists")
    else:
        df = pd.read_csv("../data/TRb_CD8_n42675_uppercase.csv", sep=";")
        # take only TCR with LD=1
        unique_tcr_epi = set(df_edge['Source']).union(set(df_edge['Target']))
        df = df[df['CDR3'].isin(unique_tcr_epi)]
        df = df[df['Epitope_peptide'].isin(unique_tcr_epi)]

        # TCR data
        # add pleiospecific information (one TCR recognizes multiple epitopes)
        df = df.set_index('CDR3')
        df_tcr_patho = pd.read_csv("../data/CDR3_Pathology.csv", sep=";", index_col='CDR3')
        df_tcr_patho[df_tcr_patho != 0] = 1 # to binary
        # if sum>1, pleiospecific
        df_tcr_patho = df_tcr_patho.sum(axis=1) > 1
        # join
        df = df.join(df_tcr_patho.to_frame(name='Pleiospecific'), on='CDR3')

        df_tcr = pd.DataFrame({'Id': df.index,
                               'Node_type': 'TCR',
                               'Pleiospecific': df['Pleiospecific']})
        df_tcr = df_tcr.drop_duplicates(subset=['Id'])
        print("Unique TCR nodes:", df_tcr.shape[0])

        # Epitope data with pathology as attribute
        df_epi = df.drop_duplicates(subset=['Epitope_peptide'])
        df_epi = pd.DataFrame({'Id': df_epi['Epitope_peptide'],
                               'Node_type': df_epi['Pathology'],
                               'Pleiospecific': 'Epitope'})
        df_epi = df_epi.drop_duplicates(subset=['Id'])
        print("Unique Epitope nodes:", df_epi.shape[0])
        # concatenate TCR-Epitope
        df = pd.concat([df_tcr, df_epi], ignore_index=True)
        print("TCR + Epitope nodes:", df.shape[0])
        df.to_csv(file, sep=';', index=False)

# # Plot correlation
# b = pd.read_csv('../data/CDR3_Pathology.csv', sep=';', index_col='CDR3')
# correlation_plot(b, 'pathology')
# b = pd.read_csv('../data/CDR3_Pathology_LD1.csv', sep=';', index_col='CDR3')
# correlation_plot(b, 'pathology_LD1')

df = create_edge_data(plot=False)
create_node_data(df)

