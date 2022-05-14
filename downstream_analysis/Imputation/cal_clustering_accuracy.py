# -*- coding: utf-8 -*-
"""
Created on Sat May 14 18:42:01 2022

@author: WBC
"""
from clustering import spectral_clustering, hierarchical_clustering,get_metric, get_metric_without_label, kmeans_clustering
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from utiles import dataNormalization
import umap



"""
####################### simulation data ################################## 
"""
datapath = "simu/simu/simu_g2000_c500_k5_s2021_"
resultpath = "simuData/simu_imputation_result/simu_g2000_c500_k5_s2021_"
dataset = ["dm1", "dm3","dm4","dm5", "dm7"]

method = ["MAGIC", "kNN-smoothing", "SAVER", "I-Impute", "SMURF_V", "SMURF_F", "SMURF_CV", "SC"]

Tool = None

for d in dataset:
    label = pd.read_csv(datapath+d+"_label.csv", header=None)
    raw = pd.read_csv(datapath+d+"_raw_count.csv", header=0, index_col=0)
    raw_label = spectral_clustering(raw.values.transpose(), k=5)
    res1 = get_metric(np.squeeze(label.values), raw_label)
    sw = get_metric_without_label(X=raw.values.transpose(), label=raw_label)
    res1["SW"] = sw['SW']
    res1 = pd.DataFrame(res1, index=["Raw"])
    res1["Dataset"] = d

    dropout = pd.read_csv(datapath+d+"_dropout_count.csv", header=0, index_col=0)
    dropout_label = spectral_clustering(dropout.values.transpose(), k=5)
    res2 = get_metric(np.squeeze(label.values), dropout_label)
    sw = get_metric_without_label(X=dropout.values.transpose(), label=dropout_label)
    res2["SW"] = sw['SW']
    res2 = pd.DataFrame(res2, index=["Dropout"])
    res2["Dataset"] = d

    res = pd.concat((res1, res2), axis=0)

    for m in method:
        impute = pd.read_csv(resultpath + d + "_"+m+".csv", header=0, index_col=0)
        impute_label = spectral_clustering(impute.values.transpose(), k=5)
        res3 = get_metric(np.squeeze(label.values), impute_label)
        sw = get_metric_without_label(X=impute.values.transpose(), label=impute_label)
        res3["SW"] = sw['SW']
        res3 = pd.DataFrame(res3, index=[m])
        res3["Dataset"] = d

        res = pd.concat((res, res3), axis=0)

    if Tool :
        ans = pd.concat((ans, res), axis=0)
    else:
        ans = res
        Tool = True


print(ans)
# ans.to_csv("simuData/simu_imputation_result/metrics.csv")


"""
####################### CellBench data ################################## 
"""
datapath = "cellbench/"
resultpath = "cellbench/"


dataset = ["p3cl", "sc_10x","sc_10x_5cl","sc_celseq2_5cl_p1","sc_celseq2_5cl_p2", "sc_celseq2_5cl_p3", "sc_celseq2", "sc_dropseq"]

clusters_n = [3,3,5,5,5,5,5,3,3]

method = ["SMURF_V", "SMURF_F", "SMURF_CV"]
method = ["MAGIC","kNN-smoothing", "SAVER","I-Impute","SMURF_V", "SMURF_F","SMURF_CV","SMURF_V_H", "SMURF_F_H","SMURF_CV_H", "SC","SC_SMURF_V", "SC_SMURF_V_H", "SC_SMURF_F", "SC_SMURF_F_H", "SC_SMURF_CV", "SC_SMURF_CV_H"]

Tool = None

for i, d in enumerate(dataset):
    label = pd.read_csv(datapath+d+"_metadata.csv", header=0, index_col=None)

    if d in ['sc_celseq2_5cl_p1', 'sc_celseq2_5cl_p2', 'sc_celseq2_5cl_p3']:
        label = label["cell_line_demuxlet"]
    else:
        label = label["cell_line"]

    raw = pd.read_csv(datapath+d+"_count_umap.csv", header=0, index_col=0)


    raw_label = spectral_clustering(raw.values, k=clusters_n[i])
    res1 = get_metric(np.squeeze(label.values), raw_label)
    sw = get_metric_without_label(X=raw.values, label=raw_label)
    res1["SW"] = sw['SW']
    res1 = pd.DataFrame(res1, index=["Raw"])
    res1["Dataset"] = d

    dropout = pd.read_csv(datapath+d+"_samp_umap.csv", header=0, index_col=0)
    # pca1 = PCA(n_components=15)
    # dropout = pca.transform(dropout.values.transpose())

    dropout_label = spectral_clustering(dropout.values, k=clusters_n[i])

    res2 = get_metric(np.squeeze(label.values), dropout_label)
    sw = get_metric_without_label(X=dropout.values, label=dropout_label)
    res2["SW"] = sw['SW']
    res2 = pd.DataFrame(res2, index=["Dropout"])
    res2["Dataset"] = d

    res = pd.concat((res1, res2), axis=0)


    for m in method:
        impute = pd.read_csv(resultpath + d + "_"+m+"_umap.csv", header=0, index_col=0)



        impute_label = spectral_clustering(impute.values, k=clusters_n[i])
        res3 = get_metric(np.squeeze(label.values), impute_label)
        sw = get_metric_without_label(X=impute.values, label=impute_label)
        res3["SW"] = sw['SW']
        res3 = pd.DataFrame(res3, index=[m])
        res3["Dataset"] = d

        res = pd.concat((res, res3), axis=0)

    if Tool :
        ans = pd.concat((ans, res), axis=0)
    else:
        ans = res
        Tool = True


print(ans)
