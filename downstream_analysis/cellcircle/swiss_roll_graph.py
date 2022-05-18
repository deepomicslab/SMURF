# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 14:23:37 2022

@author: WANG Bingchen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

"""#####################     H1-hESC     #####################"""
# methods = ["normalized","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
# # methods = [ "normalized"]
# for m in methods:
#     if m == "normalized":
#         dataFile = "test_data/H1-hESC_"+m+".csv"
#     else:
#         dataFile = "result/H1-hESC_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)

#     clusters = []
#     for col in dataFrame.columns:
#         if re.match("H1_.*", col):
#             dataFrame = dataFrame.drop([col], axis=1)
#         if re.match("G1_.*", col):
#             clusters.append(1)
#         if re.match("S_.*", col):
#             clusters.append(2)
#         if re.match("G2_.*", col):
#             clusters.append(3)


#     cellsNames = dataFrame.columns.values.tolist()
    
#     res = pd.read_csv("oval/H1-hESC_"+m+"_1D-oval.csv",header=0, index_col=0)
#     angle = res["angle(rad)"]
#     res.columns = ["angle(rad)","y","z"]
#     y = res["y"].values
#     z = res["z"].values

#     fig = plt.figure()

#     ax = plt.axes(projection="3d")

#     G1 = []
#     G2 = []
#     S = []
#     x = []
#     for i in range(len(cellsNames)):
#         if clusters[i] == 1:
#             x1 = random.gauss(0.005,0.002)
#             G1.append([x1,y[i],z[i]])
#             x.append(x1)
#         if clusters[i] == 2:
#             x2 = random.gauss(0,0.002)
#             S.append([x2,y[i],z[i]])
#             x.append(x2)
#         if clusters[i] == 3:
#             x3 = random.gauss(-0.005, 0.002)
#             G2.append([x3,y[i],z[i]])
#             x.append(x3)
#     G1 = np.array(G1)
#     G2 = np.array(G2)
#     S = np.array(S)
#     ax.scatter3D(G1[:, 0], G1[:, 1],G1[:,2])
#     ax.scatter3D(S[:, 0], S[:, 1],S[:,2])
#     ax.scatter3D(G2[:, 0], G2[:, 1],G2[:,2])
#     ax.legend(["G1", "S", "G2"])
#     plt.xlim([-0.04,0.04])
#     ax.set_title(m)
#     plt.show()
#     res.insert(1, "x", x)
#     res.to_csv("swiss_roll/H1-hESC_"+m+"_swiss_oval.csv")
    
"""###################   3Line-qPCR   #####################"""


labels = pd.read_csv("test_data/3Line-qPCR_metadata.csv", header=0, index_col=0)


methods = ["log2","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
# methods = ["SMURF_F"]
for m in methods:
    if m == "log2":
        dataFile = "test_data/3Line-qPCR_PC3_"+m+".csv"
    else:
        dataFile = "result/3Line-qPCR_PC3_"+m+".csv"
    dataFrame = pd.read_csv(dataFile, header=0, index_col=0)



    cellsNames = dataFrame.columns.values.tolist()


    res = pd.read_csv("oval/3Line-qPCR_PC3_"+m+"_1D-oval.csv",header=0, index_col=0)

    angle = res["angle(rad)"]
    res.columns = ["angle(rad)","y","z"]
    y = res["y"].values
    z = res["z"].values

    fig = plt.figure()

    ax = plt.axes(projection="3d")
    G01 = []
    G2M = []
    S = []
    x = []
    for i in range(len(angle)):
        p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
        if p == "G0/G1":
            x1 = random.gauss(0.005,0.002)
            G01.append([x1,y[i],z[i]])
            x.append(x1)
        if p == "S":
            x2 = random.gauss(0, 0.002)
            S.append([x2,y[i],z[i]])
            x.append(x2)
        if p == "G2/M":
            x3 = random.gauss(-0.005,0.002)
            G2M.append([x3,y[i],z[i]])
            x.append(x3)
    G01 = np.array(G01)
    G2M = np.array(G2M)
    S = np.array(S)
    ax.scatter3D(G01[:, 0], G01[:, 1],G01[:,2])
    ax.scatter3D(S[:, 0], S[:, 1],S[:,2])
    ax.scatter3D(G2M[:, 0], G2M[:, 1],G2M[:,2])
    ax.legend(["G0/G1", "S", "G2/M"])
    ax.set_title(m)
    plt.xlim([-0.04,0.04])
    plt.show()

    res.insert(1, "x", x)
    res.to_csv("swiss_roll/3Line-qPCR_PC3_"+m+"_swiss_oval.csv")
    
    
"""###################   mESC-Quartz   #####################"""


# labels = pd.read_csv("test_data/mESC-Quartz_metadata.csv", header=0, index_col=0)


# methods = ["FPKM_samp","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]

# for m in methods:
#     if m == "FPKM_samp":
#         dataFile = "test_data/mESC-Quartz_"+m+".csv"
#     else:
#         dataFile = "result/mESC-Quartz_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)


#     cellsNames = dataFrame.columns.values.tolist()
#     res = pd.read_csv("oval/mESC-Quartz_"+m+"_1D-oval.csv",header=0, index_col=0)
   
#     angle = res["angle(rad)"]
#     res.columns = ["angle(rad)","y","z"]
#     y = res["y"].values
#     z = res["z"].values

#     fig = plt.figure()

#     ax = plt.axes(projection="3d")


#     G01 = []
#     G2M = []
#     S = []
#     x = []
#     for i in range(len(angle)):
#         p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
#         if p == "G1":
#             x1 = random.gauss(0.005,0.002)
#             G01.append([x1,y[i],z[i]])
#             x.append(x1)
#         if p == "S":
#             x2 = random.gauss(0, 0.002)
#             S.append([x2,y[i],z[i]])
#             x.append(x2)
#         if p == "G2/M":
#             x3 = random.gauss(-0.005,0.002)
#             G2M.append([x3,y[i],z[i]])
#             x.append(x3)
#     G01 = np.array(G01)
#     G2M = np.array(G2M)
#     S = np.array(S)
#     ax.scatter3D(G01[:, 0], G01[:, 1],G01[:,2])
#     ax.scatter3D(S[:, 0], S[:, 1],S[:,2])
#     ax.scatter3D(G2M[:, 0], G2M[:, 1],G2M[:,2])
#     ax.legend(["G0/G1", "S", "G2/M"])
#     ax.set_title(m)
#     plt.xlim([-0.04,0.04])
#     plt.show()
    
#     res.insert(1, "x", x)
#     res.to_csv("swiss_roll/mESC-Quartz_"+m+"_swiss_oval.csv")

"""###################   mESC-SMARTer   #####################"""


# labels = pd.read_csv("test_data/mESC-SMARTer_metadata.csv", header=0, index_col=0)


# methods = ["count","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]

# for m in methods:
#     if m == "count":
#         dataFile = "test_data/mESC-SMARTer_"+m+".csv"
#     else:
#         dataFile = "result/mESC-SMARTer_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)


#     cellsNames = dataFrame.columns.values.tolist()
#     res = pd.read_csv("oval/mESC-SMARTer_"+m+"_1D-oval.csv",header=0, index_col=0)


#     angle = res["angle(rad)"]
#     res.columns = ["angle(rad)","y","z"]
#     y = res["y"].values
#     z = res["z"].values

#     fig = plt.figure()

#     ax = plt.axes(projection="3d")
#     G01 = []
#     G2M = []
#     S = []
#     x = []
#     for i in range(len(angle)):
#         p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
#         if p == "G1":
#             x1 = random.gauss(0.005,0.002)
#             G01.append([x1,y[i],z[i]])
#             x.append(x1)
#         if p == "S":
#             x2 = random.gauss(0, 0.002)
#             S.append([x2,y[i],z[i]])
#             x.append(x2)
#         if p == "G2M":
#             x3 = random.gauss(-0.005,0.002)
#             G2M.append([x3,y[i],z[i]])
#             x.append(x3)
#     G01 = np.array(G01)
#     G2M = np.array(G2M)
#     S = np.array(S)
#     ax.scatter3D(G01[:, 0], G01[:, 1],G01[:,2])
#     ax.scatter3D(S[:, 0], S[:, 1],S[:,2])
#     ax.scatter3D(G2M[:, 0], G2M[:, 1],G2M[:,2])
#     ax.legend(["G0/G1", "S", "G2/M"])
#     ax.set_title(m)
#     plt.xlim([-0.04,0.04])
#     ax.set_title(m)
#     plt.show()
#     res.insert(1, "x", x)
#     res.to_csv("swiss_roll/mESC-SMARTer_"+m+"_swiss_oval.csv")
# #
