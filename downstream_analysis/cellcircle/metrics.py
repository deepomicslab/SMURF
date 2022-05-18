# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:00:27 2022

@author: WANG Bingchen
"""


import pandas as pd
import numpy as np
from clustering import spectral_clustering, kmeans_clustering, hierarchical_clustering
from sklearn import metrics


def calc_change_index(name, angle, label, n_cluster=3):
    N = len(name)
    comb = zip(name, angle, label)
    order = sorted(comb, key=lambda comb: comb[1])
    orderName, orderAngle, orderLabel = list(zip(*order))[0],list(zip(*order))[1],list(zip(*order))[2]
    s_c = 0
    label_pre = orderLabel[0]
    for i in range(1,N):
        l = orderLabel[i]
        if l != label_pre:
            s_c += 1
            label_pre = l
    if orderLabel[0] != orderLabel[-1]:
        s_c += 1
    score = 1 - ((s_c - n_cluster + 1)/(N - n_cluster))
    return score


def calc_pcc(name, pos, label, n_cluster=3):
    newlabels = kmeans_clustering(pos.transpose(), k=n_cluster)
    res = []
    for i in range(n_cluster):
        newlabels += 1
        newlabels[newlabels == n_cluster+1] = 1
        res.append(np.min(np.corrcoef(newlabels, label)))

        newlabels_inv = -(newlabels - (n_cluster + 1))
        res.append(np.min(np.corrcoef(newlabels_inv, label)))
    return np.max(res)









if __name__ == '__main__':

    data = pd.read_csv("result/ovalDataTest.csv", header=0, index_col=0)
    data = pd.read_csv('test_data/GSE64016_H1andFUCCI_normalized_EC_no_imputation_cell_circle.csv', header=0, index_col=0)
    data = pd.read_csv('test_data/H1-hESC_normalized_1D-oval.csv', header=0,
                       index_col=0)
    cellsName = np.array(data.columns)
    cellsAngle = data.values[0, :]
    cellsPos = data.values[1:3, :]

    def buildLabel(str):
        if str == "G1":
            return 1
        elif str == "G2":
            return 2
        else:
            return 3

    cellsLabel = np.array([buildLabel(x[:2]) for x in cellsName])
    a = calc_change_index(cellsName, cellsAngle, cellsLabel)
    b = calc_pcc(cellsName, cellsPos, cellsLabel)
    print(a)
    print(b)