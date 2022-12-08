import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from sklearn import decomposition
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans



if __name__ == "__main__":
    
    geno_file = sys.argv[1]
    coord_file = sys.argv[2]

    ### read data
    data = pd.read_csv(geno_file, low_memory=False).dropna(axis=0).drop(["CHROM", "POS", "REF", "ALT"], axis=1)

    ### read coord file
    coord_df  = pd.read_csv(coord_file, sep="\t", lineterminator='\n')
    coord_ind = list(coord_df.loc[:, "UUID"])

    ### extract coordinated data
    named_data = data.loc[:, data.columns.intersection(coord_ind)].transpose()
    named_coord = coord_df.loc[coord_df["UUID"].isin(named_data.index), :]
    

    ### add unamed data
    unnamed_data = data.drop(named_data.index, axis=1).transpose()
    # unnamed_classes = [5] * len(unnamed_data.index)
    unnamed_coord = pd.DataFrame(data={"UUID":unnamed_data.index,
                                       "Superpopulation code":["Unknown"]*len(unnamed_data.index)})
    ### complete data
    complete_data = pd.concat([named_data, unnamed_data], axis=0)
    coord_df = pd.concat([named_coord, unnamed_coord])
    coord_df = coord_df.set_index("UUID").reindex(index=complete_data.index)
    
    coord_classes, numeric_classes = np.unique(coord_df["Superpopulation code"], return_inverse=True)
    coord_df["Superpopulation numeric"] = numeric_classes

    ### compute pca
    pca = decomposition.PCA(n_components=10)
    pca.fit(complete_data)
    pca_val = pca.transform(complete_data)
    ##################################################

    ### compute hierarchical clustering
    hierarchical_cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
    coord_df["hierarchical"] = hierarchical_cluster.fit_predict(complete_data)
    ##################################################

    ### compute k means
    km = KMeans(
        n_clusters=5, init='random',
        n_init=10, max_iter=300, 
        tol=1e-04, random_state=0
        )
    coord_df["kmeans"] = km.fit_predict(complete_data)
    #############################################

    ### plot pca custom
    cols = ["gray",  "#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"]

    fig = plt.figure(1, figsize=(12, 9))
    
    ### pca plot
    sub_pca = fig.add_subplot(1, 3, 1)
    sub_pca.title.set_text('PCA vals')
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "Unknown", 0], 
                    pca_val[coord_df["Superpopulation code"] == "Unknown", 1], color=cols[0], edgecolor=None, alpha=0.5, s=15)
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "AFR", 0], 
                    pca_val[coord_df["Superpopulation code"] == "AFR", 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "AMR", 0], 
                    pca_val[coord_df["Superpopulation code"] == "AMR", 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "EAS", 0], 
                    pca_val[coord_df["Superpopulation code"] == "EAS", 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "EUR", 0], 
                    pca_val[coord_df["Superpopulation code"] == "EUR", 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
    sub_pca.scatter(pca_val[coord_df["Superpopulation code"] == "SAS", 0], 
                    pca_val[coord_df["Superpopulation code"] == "SAS", 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
    sub_pca.legend(coord_df["Superpopulation code"])

    ### hierarchical plot
    sub_hier = fig.add_subplot(1, 3, 2)
    sub_hier.title.set_text("Hierarchical clustering")
    sub_hier.scatter(pca_val[coord_df["hierarchical"] == 0, 0],
                     pca_val[coord_df["hierarchical"] == 0, 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
    sub_hier.scatter(pca_val[coord_df["hierarchical"] == 1, 0],
                     pca_val[coord_df["hierarchical"] == 1, 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
    sub_hier.scatter(pca_val[coord_df["hierarchical"] == 2, 0],
                     pca_val[coord_df["hierarchical"] == 2, 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
    sub_hier.scatter(pca_val[coord_df["hierarchical"] == 3, 0],
                     pca_val[coord_df["hierarchical"] == 3, 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
    sub_hier.scatter(pca_val[coord_df["hierarchical"] == 4, 0],
                     pca_val[coord_df["hierarchical"] == 4, 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
    # sub_hier.legend(coord_df["Superpopulation code"])

    ### kmeans plot
    sub_kms = fig.add_subplot(1, 3, 3)
    sub_kms.title.set_text("KMeans clustering")
    sub_kms.scatter(pca_val[coord_df["kmeans"] == 0, 0],
                    pca_val[coord_df["kmeans"] == 0, 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
    sub_kms.scatter(pca_val[coord_df["kmeans"] == 1, 0],
                    pca_val[coord_df["kmeans"] == 1, 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
    sub_kms.scatter(pca_val[coord_df["kmeans"] == 2, 0],
                    pca_val[coord_df["kmeans"] == 2, 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
    sub_kms.scatter(pca_val[coord_df["kmeans"] == 3, 0],
                    pca_val[coord_df["kmeans"] == 3, 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
    sub_kms.scatter(pca_val[coord_df["kmeans"] == 4, 0],
                    pca_val[coord_df["kmeans"] == 4, 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
    # sub_kms.legend(coord_df["Superpopulation code"])
                     
    plt.show()
