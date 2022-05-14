import pandas as pd

import iimpute

data = pd.read_csv("traj/samp/cellmix4.csv", index_col=0, header=0)

# create I-Impute object
# iimpute_operator = iimpute.IImpute(normalize=True, iteration=False)
#
# # impute
# imputed_data = iimpute_operator.impute(data)
#
# # store result to a file
# imputed_data.to_csv('cellbench/sc_10x_5cl_C-Impute.csv')

# iterative mode
iimpute_operator = iimpute.IImpute(iteration=True, normalize=True, n=5)

# impute
imputed_data = iimpute_operator.impute(data)

# store result to a file
imputed_data.to_csv('traj/I-Impute/cellmix4.csv')



