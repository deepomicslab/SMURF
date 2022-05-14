# -*- coding: utf-8 -*-
"""
Created on Sat May 14 21:44:00 2022

@author: WBC
"""

import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from scipy import stats
from math import sqrt



def calculateRMSE(data1, data2):

    return sqrt(mean_squared_error(data1, data2))


datapath = "simu/simu/simu_g2000_c500_k5_s2021_"
resultpath = "simuData/simu_imputation_result/simu_g2000_c500_k5_s2021_"
dataSets = ["dm1", "dm3","dm4","dm5", "dm7"]

methods = ["MAGIC", "kNN-smoothing",  "SAVER", "I-Impute","SMURF_V", "SMURF_F", "SMURF_CV"]


imputation_rmse = np.zeros((len(methods), len(dataSets)))

for j, ds in enumerate(dataSets):
    origin_data = pd.read_csv(datapath+ds + '_raw_count.csv', index_col=0, header=0)
    for i, m in enumerate(methods):
        impute_data = pd.read_csv(resultpath + ds +m+'.csv', index_col=0, header=0)

        imputation_rmse[i, j] = calculateRMSE(origin_data.values, impute_data.values)


imputation_rmse_dataframe = pd.DataFrame(imputation_rmse, index=methods, columns=dataSets)
print(imputation_rmse_dataframe)




imputation = np.zeros((len(methods), len(dataSets)))

for j, ds in enumerate(dataSets):
    origin_data = pd.read_csv(datapath+ds + '_raw_count.csv', index_col=0, header=0).transpose()


    for i, m in enumerate(methods):
        impute_data = pd.read_csv(resultpath + ds +m+'.csv', index_col=0, header=0).transpose()


        rows, cols = np.shape(origin_data)
        a = []

        for x in origin_data.columns:
            corr = stats.pearsonr(origin_data[x], impute_data[x])[0]
            a.append(corr)
        print(a)

        # imputation[i, j] = [np.median(a), np.quantile(a, 0.25), np.quantile(a, 0.75)]
        imputation[i, j] = np.quantile(a, 0.25)



correlation = pd.DataFrame(imputation, index=methods, columns=dataSets)

print(correlation)