# -*- coding: utf-8 -*-
"""
Created on Sat May 14 17:19:13 2022

@author: WANG Bingchen
"""

import pandas as pd
import numpy as np
from metrics import calc_pcc, calc_change_index

datasets = ["H1-hESC", "3Line-qPCR", "mESC-Quartz", "mESC-SMARTer"]
methods = ["KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F","SMURF_F_H"]

"""###############   H1-hESC   ################"""
methods = ["normalized","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F","SMURF_F_H"]
#
dic = {}
for m in methods:
    data = pd.read_csv("oval/H1-hESC_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/H1-hESC_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel1(str):
        if str == "G1":
            return 1
        elif str == "G2":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel1(x[:2]) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "H1-hESC"]
res1 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()


"""###############   3Line-qPCR   ################"""
methods = ["log2","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F","SMURF_F_H"]


dic = {}
labels = pd.read_csv("test_data/3Line-qPCR_metadata.csv", header=0, index_col=0)
for m in methods:
    data = pd.read_csv("oval/3Line-qPCR_H9_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/3Line-qPCR_H9_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel2(str):
        p = labels[labels["Cell"] == str]["Phase"].values[0]
        if p == "G0/G1":
            return 1
        elif p == "G2/M":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel2(x) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "3Line-qPCR_H9"]
res2 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()



dic = {}
labels = pd.read_csv("test_data/3Line-qPCR_metadata.csv", header=0, index_col=0)
for m in methods:
    data = pd.read_csv("oval/3Line-qPCR_MB_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/3Line-qPCR_MB_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel2(str):
        p = labels[labels["Cell"] == str]["Phase"].values[0]
        if p == "G0/G1":
            return 1
        elif p == "G2/M":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel2(x) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "3Line-qPCR_MB"]
res5 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()
print(res5)



dic = {}
labels = pd.read_csv("test_data/3Line-qPCR_metadata.csv", header=0, index_col=0)
for m in methods:
    data = pd.read_csv("oval/3Line-qPCR_PC3_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/3Line-qPCR_PC3_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel2(str):
        p = labels[labels["Cell"] == str]["Phase"].values[0]
        if p == "G0/G1":
            return 1
        elif p == "G2/M":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel2(x) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "3Line-qPCR_PC3"]
res6 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()
print(res6)


"""###############   mESC-Quartz   ################"""
methods = ["FPKM_samp","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F","SMURF_F_H"]

dic = {}
labels = pd.read_csv("test_data/mESC-Quartz_metadata.csv", header=0, index_col=0)
for m in methods:
    data = pd.read_csv("oval/mESC-Quartz_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/mESC-Quartz_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel2(str):
        p = labels[labels["Cell"] == str]["Phase"].values[0]
        if p == "G1":
            return 1
        elif p == "G2/M":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel2(x) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "mESC-Quartz"]
res3 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()


"""###############   mESC-SMARTer   ################"""
methods = ["count","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F","SMURF_F_H"]

dic = {}
labels = pd.read_csv("test_data/mESC-SMARTer_metadata.csv", header=0, index_col=0)
for m in methods:
    data = pd.read_csv("oval/mESC-SMARTer_"+m+"_1D-oval.csv", header=0, index_col=0)
    data = data.transpose()
    # data.to_csv("oval/mESC-SMARTer_"+m+"_1D-oval.csv")
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]


    def buildLabel2(str):
        p = labels[labels["Cell"] == str]["Phase"].values[0]
        if p == "G1":
            return 1
        elif p == "G2/M":
            return 2
        else:
            return 3


    cellsLabel = np.array([buildLabel2(x) for x in cellsName])
    CI = calc_change_index(cellsName, cellsAngle, cellsLabel)
    PCC = calc_pcc(cellsName, cellsPos, cellsLabel)

    dic[m] = [CI, PCC, "mESC-SMARTer"]
res4 = pd.DataFrame(dic, index=["Change Index", "PCC", "Data Set"]).transpose()


res = pd.concat((res1, res2,res5,res6, res3, res4), axis=0)
res.to_csv("oval/metrics.csv")
