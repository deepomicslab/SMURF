import magic
import pandas as pd



# dataSets = ["sc_celseq2", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2", "sc_celseq2_5cl_p3", "sc_dropseq"]

X = pd.read_csv("traj/samp/RNAmix_celseq2.csv", index_col=0, header=0)
X = X.transpose()
print(X.head())
magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(X, genes="all_genes", t_max=30)
X_magic = X_magic.transpose()
X_magic.to_csv('traj/MAGIC/RNAmix_celseq2.csv')


# dataSets = ["dm1", "dm3", "dm4","dm5", "dm7"]
# for ds in dataSets:
#     origin_data = pd.read_csv("simu/simu/simu_g2000_c500_k5_s2021_" + ds + "_dropout_count.csv", index_col=0, header=0)
#     origin_data = origin_data.transpose()
#
#     magic_operator = magic.MAGIC()
#     impute_data = magic_operator.fit_transform(origin_data)
#     impute_data = impute_data.transpose()
#
#     impute_data.to_csv("simuData/simu_imputation_result/simu_g2000_c500_k5_s2021_" + ds + "_MAGIC.csv")






#
# for ds in dataSets:
#     origin_data = pd.read_csv("cellbench/" + ds + "_samp.csv", index_col=0, header=0)
#     origin_data = origin_data.transpose()
#
#     magic_operator = magic.MAGIC()
#     impute_data = magic_operator.fit_transform(origin_data)
#     impute_data = impute_data.transpose()
#
#     impute_data.to_csv("cellbench_result/" + ds + "_MAGIC.csv")
