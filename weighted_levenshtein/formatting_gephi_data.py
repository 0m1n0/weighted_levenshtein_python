import pandas as pd
from pathlib import Path
pd.options.mode.chained_assignment = None  # default='warn'


# formatting data for Gephi network visualization software
def rescale_wld_to_similarity(x):
    # for LD=1, we have WLD between 8-15
    wld_list = list(range(8, 16))
    # rescaled: between 1-8 as similarity score
    rescaled_value = x - 2 * wld_list.index(x)
    return rescaled_value


def create_edge_data():
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
        # TCR-Epitope edges with fixed max weight
        df_epi = pd.read_csv("../data/TRb_CD8_n42675_uppercase.csv", sep=";")
        # keep only TCR with LD=1
        df_epi = df_epi[df_epi['CDR3'].isin(tcr_list)]
        df_epi = pd.DataFrame({'Source': df_epi['CDR3'],
                               'Target': df_epi['Epitope_peptide'],
                               'Weight': 8,
                               'Edge_type': 'TCR-Epitope'})
        # concatenate TCR-Epitope
        df = pd.concat([df_tcr, df_epi], ignore_index=True)
        print("TCR + Epitope edges:", df.shape[0])
        df.to_csv(file, sep=';', index=False)
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
        df_tcr = pd.DataFrame({'Id': df['CDR3'],
                               'Node_type': 'TCR'})
        df_tcr = df_tcr.drop_duplicates(subset=['Id'])
        print("Unique TCR nodes:", df_tcr.shape[0])
        # Epitope data with pathology as attribute
        df_epi = df.drop_duplicates(subset=['Epitope_peptide'])
        df_epi = pd.DataFrame({'Id': df_epi['Epitope_peptide'],
                               'Node_type': df_epi['Pathology']})
        df_epi = df_epi.drop_duplicates(subset=['Id'])
        print("Unique Epitope nodes:", df_epi.shape[0])
        # concatenate TCR-Epitope
        df = pd.concat([df_tcr, df_epi], ignore_index=True)
        print("TCR + Epitope nodes:", df.shape[0])
        df.to_csv(file, sep=';', index=False)


df = create_edge_data()
create_node_data(df)
