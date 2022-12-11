# _*_ coding: utf-8 _*_
"""
Time:     2021/11/29 12:19
Author:   WANG Bingchen
Version:  V 0.1
File:     selfConsistency.py
Describe: 
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utiles import dataNormalization



def initialSCMatrix(cell_num):
    weight0 = float(1)/float(cell_num-1)
    crm_value = np.log(weight0/(1-weight0))
    cell_relation_Matrix = np.ones((cell_num, cell_num))*crm_value
    # randWeight = np.random.random([cell_num, cell_num])
    # cell_relation_Matrix = cell_relation_Matrix + crm_value*randWeight

    row, col = np.diag_indices_from(cell_relation_Matrix)
    cell_relation_Matrix[row, col] = float("-inf")

    return cell_relation_Matrix






# def SC_weight(cell_relation_Matrix):
#     cell_relation_Matrix = np.exp(cell_relation_Matrix)
#     return cell_relation_Matrix

def sigmoid(x):
    x = 1/(1+np.exp(-1*x))
    return x

def d_sigmoid(r):
    # return sigmoid(x)*(1-sigmoid(x))
    return r*(1-r)


def objective_fun(n, A, r):
    size = r.shape[0]
    fun = 0.1*np.dot(r.transpose(), r).trace()/(size**2)
    f1 = A
    f2 = A
    for i in range(n):
        f1 = np.dot(r, f1)
        if i == 0:
            f = f1 - f2
        else:
            f2 = np.dot(r, f2)
            f = (0.3)*(f1 - f2)
        fun += np.dot(f.transpose(), f).trace()/(size**2)
    return 0.5*fun


def grad_objective_fun(n, A, r):
    size = r.shape[0]
    grad_f = 0.1*r*d_sigmoid(r)/(size**2)
    grad_f1 = A
    grad_f2 = A

    grad_f3 = np.dot(A.transpose(), (r - np.identity(r.shape[0])))
    for i in range(n):
        grad_f1 = np.dot(r, grad_f1)
        if i == 0:
            d = grad_f1 - grad_f2
            f = np.dot(d, A.transpose()) * d_sigmoid(r)

        else:
            grad_f2 = np.dot(r, grad_f2)
            d = grad_f1 - grad_f2
            f = (0.3)*np.dot(d, grad_f3) * d_sigmoid(r)
            grad_f3 = np.dot(grad_f3, r)

        grad_f += f/(size**2)
    return grad_f


def sc_optimize(ln_r, data, steps, alpha, eps, n):
    losses = []
    loss_pre = 1e20
    for i in range(steps):
        r = sigmoid(ln_r)

        loss = objective_fun(n, data, r)
        losses.append(loss)
        if np.abs(loss - loss_pre) < eps:
            break
        else:
            grad = grad_objective_fun(n, data, r)
            ln_r = ln_r - alpha*grad

            ln_r[ln_r >= 0] = -(np.e)/100
            # ln_r = (ln_r + ln_r.transpose())/2
        if (i+1) % 10 == 0:
            print("step", i+1)
    return ln_r, sigmoid(ln_r), np.dot(r, data), losses



# dataFrame = pd.read_csv("cellbench/sc_10x_samp.csv", header=0, index_col=0).astype("float64")
# crm = initialSCMatrix(dataFrame.shape[1])
#
# ln_r, r, res, l = sc_optimize(crm, dataFrame.values.transpose(), 300, 1e-6, 1, 10)
# plt.plot(l)
# plt.show()
# print(dataFrame.values.transpose())
# print(res)

if __name__ == "__main__":
    dataFrame = pd.read_csv("traj/RNAmix_sortseq/genebycell.csv", header=0, index_col=0)
    dataFrame = pd.read_csv("traj/samp/RNAmix_sortseq.csv", header=0, index_col=0)
    print(dataFrame)
    df, sizefactor = dataNormalization(dataFrame)
    crm = initialSCMatrix(dataFrame.shape[1])



    ln_r, r, res, l = sc_optimize(crm, df.transpose().values, 3000,0.1, 1, 3)
    plt.plot(l)
    plt.show()

    newDataFrame = pd.DataFrame(res.transpose(), index=dataFrame.index, columns=dataFrame.columns)
    newDataFrame = newDataFrame * sizefactor
    print(newDataFrame)

    newDataFrame.to_csv("traj/SC/RNAmix_sortseq.csv")
    # newDataFrame.to_csv("simuData/simu_imputation_result/simu_g2000_c500_k5_s2021_dm1_SC.csv")













