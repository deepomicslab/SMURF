from knn_smooth import knn_smoothing
import pandas as pd
import numpy as np

origin_data = pd.read_csv("traj/samp/RNAmix_sortseq.csv", index_col=0, header=0).astype(np.float64)
print(origin_data)
impute_data = knn_smoothing(origin_data.values, k=5)

impute_dataframe = pd.DataFrame(impute_data, index=origin_data.index, columns=origin_data.columns)
impute_dataframe.to_csv('traj/kNN-smoothing/RNAmix_sortseq.csv')
print(impute_dataframe)
# dataSets = ["sc_celseq2", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2", "sc_celseq2_5cl_p3", "sc_dropseq"]


# for ds in dataSets:
#     origin_data = pd.read_csv("cellbench/" + ds + "_samp.csv", index_col=0, header=0).astype(np.float64)
#     # p, n = origin_data.shape
#     # print('The expression matrix contains %d genes and %d cells.' % (p, n))
#     impute_data = knn_smoothing(origin_data.values, k=5, d=10)
#
#     impute_dataframe = pd.DataFrame(impute_data, index=origin_data.index, columns=origin_data.columns)
# #
#     impute_dataframe.to_csv("cellbench_result/" + ds + "_knnsmooth.csv")


#
# dataSets = ["dm1", "dm3", 'dm4',"dm5", "dm7"]
# for ds in dataSets:
#     origin_data = pd.read_csv("simu/simu/simu_g2000_c500_k5_s2021_" + ds + "_dropout_count.csv", index_col=0, header=0).astype(np.float64)
#     impute_data = knn_smoothing(origin_data.values,k=3)
#
#     impute_dataframe = pd.DataFrame(impute_data, index=origin_data.index, columns=origin_data.columns)
#     impute_dataframe.to_csv("simuData/simu_imputation_result/simu_g2000_c500_k5_s2021_" + ds + "_kNN-smoothing.csv")
