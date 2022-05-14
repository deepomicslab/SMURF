# -*- coding: utf-8 -*-
"""
Created on Sat May 14 17:38:20 2022

@author: WANG Bingchen
"""

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from cell_circle import CellCircle
import smurf


"""#####################     H1-hESC     #####################"""
# methods = ["KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
# methods = [ "normalized"]
# for m in methods:
#     dataFile = "result/H1-hESC_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)
#
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
#
#
#     cellsNames = dataFrame.columns.values.tolist()
#     cell_circle_mapper = CellCircle()
#
#     res = cell_circle_mapper.cal_cell_circle(dataFrame.values, a=3, b=2, k=0.2)
#
#     angle = res["angle"]
#
#     plane_embedding = res["plane_embedding"]
#
#     fig = plt.figure(1)
#     fig.clf
#     ax = fig.add_subplot(111)
#     G1 = []
#     G2 = []
#     S = []
#     for i in range(len(plane_embedding)):
#         if clusters[i] == 1:
#             G1.append(plane_embedding[i, :])
#         if clusters[i] == 2:
#             S.append(plane_embedding[i, :])
#         if clusters[i] == 3:
#             G2.append(plane_embedding[i, :])
#     G1 = np.array(G1)
#     G2 = np.array(G2)
#     S = np.array(S)
#     ax.scatter(G1[:, 0], G1[:, 1], s=1)
#     ax.scatter(S[:, 0], S[:, 1], s=1)
#     ax.scatter(G2[:, 0], G2[:, 1], s=1)
#     ax.legend(["G1", "S", "G2"])
#     ax.set_title(m)
#     plt.show()
#
#     finalcoordinate = pd.DataFrame(np.array(angle).transpose(), columns=cellsNames, index=["angle(rad)"])
#     finalcoordinate.loc["x"] = np.squeeze(plane_embedding.transpose()[0, :])
#     finalcoordinate.loc["y"] = np.squeeze(plane_embedding.transpose()[1, :])
#     finalcoordinate.to_csv("oval/H1-hESC_"+m+"_1D-oval.csv")


"""###################   3Line-qPCR   #####################"""


labels = pd.read_csv("test_data/3Line-qPCR_metadata.csv", header=0, index_col=0)


methods = ["log2","KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
methods = ["SMURF_F"]
for m in methods:
    if m == "log2":
        dataFile = "test_data/3Line-qPCR_MB_"+m+".csv"
    else:
        dataFile = "result/3Line-qPCR_MB_"+m+".csv"
    dataFrame = pd.read_csv(dataFile, header=0, index_col=0)



    cellsNames = dataFrame.columns.values.tolist()
    cell_circle_mapper = CellCircle()

    res = cell_circle_mapper.cal_cell_circle(dataFrame.values, a=3, b=2, k=0.2)

    angle = res["angle"]

    plane_embedding = res["plane_embedding"]

    fig = plt.figure(1)
    fig.clf
    ax = fig.add_subplot(111)
    G01 = []
    G2M = []
    S = []
    for i in range(len(plane_embedding)):
        p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
        if p == "G0/G1":
            G01.append(plane_embedding[i, :])
        if p == "S":
            S.append(plane_embedding[i, :])
        if p == "G2/M":
            G2M.append(plane_embedding[i, :])
    G01 = np.array(G01)
    G2M = np.array(G2M)
    S = np.array(S)
    ax.scatter(G01[:, 0], G01[:, 1], s=1)
    ax.scatter(S[:, 0], S[:, 1], s=1)
    ax.scatter(G2M[:, 0], G2M[:, 1], s=1)
    ax.legend(["G0/G1", "S", "G2/M"])
    ax.set_title(m)
    plt.show()

    finalcoordinate = pd.DataFrame(np.array(angle).transpose(), columns=cellsNames, index=["angle(rad)"])
    finalcoordinate.loc["x"] = np.squeeze(plane_embedding.transpose()[0, :])
    finalcoordinate.loc["y"] = np.squeeze(plane_embedding.transpose()[1, :])
    finalcoordinate.transpose()
    finalcoordinate.to_csv("oval/3Line-qPCR_MB_"+m+"_1D-oval.csv")


"""###################   mESC-Quartz   #####################"""

#
# labels = pd.read_csv("test_data/mESC-Quartz_metadata.csv", header=0, index_col=0)
#

# methods = ["KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
# methods = ["FPKM_samp"]
# for m in methods:
#     if m == "FPKM_samp":
#         dataFile = "test_data/mESC-Quartz_"+m+".csv"
#     else:
#         dataFile = "result/mESC-Quartz_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)
#
#
#     cellsNames = dataFrame.columns.values.tolist()
#     cell_circle_mapper = CellCircle(n_neighbors=4)
#
#     res = cell_circle_mapper.cal_cell_circle(dataFrame.values, a=3, b=2, k=0.2)
#
#     angle = res["angle"]
#
#     plane_embedding = res["plane_embedding"]
#
#     fig = plt.figure(1)
#     fig.clf
#     ax = fig.add_subplot(111)
#     G01 = []
#     G2M = []
#     S = []
#     for i in range(len(plane_embedding)):
#         p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
#         if p == "G1":
#             G01.append(plane_embedding[i, :])
#         if p == "S":
#             S.append(plane_embedding[i, :])
#         if p == "G2/M":
#             G2M.append(plane_embedding[i, :])
#     G01 = np.array(G01)
#     G2M = np.array(G2M)
#     S = np.array(S)
#     ax.scatter(G01[:, 0], G01[:, 1], s=10)
#     ax.scatter(S[:, 0], S[:, 1], s=10)
#     ax.scatter(G2M[:, 0], G2M[:, 1], s=10)
#     ax.legend(["G1", "S", "G2/M"])
#     ax.set_title(m)
#     plt.show()
#
#     finalcoordinate = pd.DataFrame(np.array(angle).transpose(), columns=cellsNames, index=["angle(rad)"])
#     finalcoordinate.loc["x"] = np.squeeze(plane_embedding.transpose()[0, :])
#     finalcoordinate.loc["y"] = np.squeeze(plane_embedding.transpose()[1, :])
#     finalcoordinate.to_csv("oval/mESC-Quartz_"+m+"_1D-oval.csv")


"""###################   mESC-SMARTer   #####################"""


# labels = pd.read_csv("test_data/mESC-SMARTer_metadata.csv", header=0, index_col=0)
#
#
# methods = ["KNN-smoothing", "MAGIC", "SAVER","C-Impute" ,"I-Impute", "SMURF_F", "SMURF_F_H"]
# methods = [ "count"]
# for m in methods:
#     if m == "count":
#         dataFile = "test_data/mESC-SMARTer_"+m+".csv"
#     else:
#         dataFile = "result/mESC-SMARTer_"+m+".csv"
#     dataFrame = pd.read_csv(dataFile, header=0, index_col=0)
#
#
#     cellsNames = dataFrame.columns.values.tolist()
#     cell_circle_mapper = CellCircle()
#
#     res = cell_circle_mapper.cal_cell_circle(dataFrame.values, a=3, b=2, k=0.2)
#
#     angle = res["angle"]
#
#     plane_embedding = res["plane_embedding"]
#
#     fig = plt.figure(1)
#     fig.clf
#     ax = fig.add_subplot(111)
#     G01 = []
#     G2M = []
#     S = []
#     for i in range(len(plane_embedding)):
#         p = labels[labels["Cell"] == cellsNames[i]]["Phase"].values[0]
#         if p == "G1":
#             G01.append(plane_embedding[i, :])
#         if p == "S":
#             S.append(plane_embedding[i, :])
#         if p == "G2M":
#             G2M.append(plane_embedding[i, :])
#     G01 = np.array(G01)
#     G2M = np.array(G2M)
#     S = np.array(S)
#     ax.scatter(G01[:, 0], G01[:, 1], s=1)
#     ax.scatter(S[:, 0], S[:, 1], s=1)
#     ax.scatter(G2M[:, 0], G2M[:, 1], s=1)
#     ax.legend(["G1", "S", "G2/M"])
#     ax.set_title(m)
#     plt.show()
#
#     finalcoordinate = pd.DataFrame(np.array(angle).transpose(), columns=cellsNames, index=["angle(rad)"])
#     finalcoordinate.loc["x"] = np.squeeze(plane_embedding.transpose()[0, :])
#     finalcoordinate.loc["y"] = np.squeeze(plane_embedding.transpose()[1, :])
#     finalcoordinate.to_csv("oval/mESC-SMARTer_"+m+"_1D-oval.csv")