# _*_ coding: utf-8 _*_
"""
Time:     2021/7/22 17:34
Author:   WANG Bingchen
Version:  V 0.1
File:     optimize.py
Describe: 
"""


from scipy.special import psi
from numpy import log


def CVOptimize(A, P, Q, mu, Vg, Vc, v, g, c):
    DFP1 = (2 * mu[g][c] / v[g][c]) * Q[:, c] * log(mu[g][c]) + mu[g][c] * Q[:, c] / v[g][c]
    DFP2 = -2 * mu[g][c] * log(v[g][c]) * Q[:, c] / v[g][c]
    DFP3 = -psi((mu[g][c] ** 2) / v[g][c]) * 2 * mu[g][c] * Q[:, c] / v[g][c]
    DFP4 = psi(A[g][c] + ((mu[g][c] ** 2) / v[g][c])) * 2 * mu[g][c] * Q[:, c] / v[g][c]
    DFP5 = -(2 * mu[g][c] * Q[:, c] / v[g][c]) * log(1 + mu[g][c] / v[g][c]) - (
            A[g][c] + mu[g][c] ** 2 / v[g][c]) * Q[:, c] / (v[g][c] + mu[g][c])
    DFP = DFP1 + DFP2 + DFP3 + DFP4 + DFP5

    DFQ1 = (2 * mu[g][c] / v[g][c]) * P[g, :] * log(mu[g][c]) + mu[g][c] * P[g, :] / v[g][c]
    DFQ2 = -2 * mu[g][c] * log(v[g][c]) * P[g, :] / v[g][c]
    DFQ3 = -psi((mu[g][c] ** 2) / v[g][c]) * 2 * mu[g][c] * P[g, :] / v[g][c]
    DFQ4 = psi(A[g][c] + mu[g][c] ** 2 / v[g][c]) * 2 * mu[g][c] * P[g, :] / v[g][c]
    DFQ5 = -(2 * mu[g][c] * P[g, :] / v[g][c]) * log(1 + mu[g][c] / v[g][c]) - (
            A[g][c] + mu[g][c] ** 2 / v[g][c]) * P[g, :] / (v[g][c] + mu[g][c])
    DFQ = DFQ1 + DFQ2 + DFQ3 + DFQ4 + DFQ5

    DFVg1 = -(mu[g][c] ** 2) * log(mu[g][c]) / Vg[g][0] ** 2
    DFVg2 = (mu[g][c] ** 2) * log(Vg[g][0]) / (Vg[g][0] ** 2) - (mu[g][c] ** 2) / (Vg[g][0] ** 2)
    DFVg3 = psi((mu[g][c] ** 2) / Vg[g][0]) * ((mu[g][c] ** 2) / Vg[g][0] ** 2)
    DFVg4 = -psi(A[g][c] + mu[g][c] ** 2 / Vg[g][0]) * ((mu[g][c] ** 2) / Vg[g][0] ** 2)
    DFVg5 = ((mu[g][c] ** 2) / (Vg[g][0] ** 2)) * log(1 + mu[g][c] / Vg[g][0]) + (
            A[g][c] + mu[g][c] ** 2 / Vg[g][0]) * (mu[g][c] / (Vg[g][0] ** 2 + mu[g][c] * Vg[g][0]))
    DFVg = DFVg1 + DFVg2 + DFVg3 + DFVg4 + DFVg5

    return DFP, DFQ, DFVg


def FanoOptimize(A, P, Q, mu, bg, bc, b, g, c):
    DFP1 = -Q[:, c] / b[g][c] * log(b[g][c])
    DFP2 = -psi(mu[g][c] / b[g][c]) * Q[:, c] / b[g][c]
    DFP3 = psi(A[g][c] + mu[g][c] / b[g][c]) * Q[:, c] / b[g][c]
    DFP4 = -(1 / b[g][c]) * log(1 + 1 / b[g][c]) * Q[:, c]
    DFP = DFP1 + DFP2 + DFP3 + DFP4

    DFQ1 = -P[g, :] / b[g][c] * log(b[g][c])
    DFQ2 = -psi(mu[g][c] / b[g][c]) * P[g][:] / b[g][c]
    DFQ3 = psi(A[g][c] + mu[g][c] / b[g][c]) * P[g][:] / b[g][c]
    DFQ4 = -(1 / b[g][c]) * log(1 + 1 / b[g][c]) * P[g, :]
    DFQ = DFQ1 + DFQ2 + DFQ3 + DFQ4

    DFBg1 = mu[g][c] / (bg[g][0] ** 2) * (log(bg[g][0]) - 1)
    DFBg2 = psi(mu[g][c] / bg[g][0]) * mu[g][c] / (bg[g][0] ** 2)
    DFBg3 = -psi(A[g][c] + mu[g][c] / bg[g][0]) * mu[g][c] / (bg[g][0] ** 2)
    DFBg4 = -(mu[g][c] * log(1 + 1 / bg[g][0]) / (bg[g][0] ** 2) + (A[g][c] * bg[g][0] + mu[g][c]) / (
            (1 + bg[g][0]) * (bg[g][0] ** 2)))
    DFBg = DFBg1 + DFBg2 + DFBg3 + DFBg4

    return DFP, DFQ, DFBg


def CCVOptimize(A, P, Q, mu, ag, ac, a, g, c):
    DFP1 = -Q[:, c] / (a[g][c] * mu[g][c])
    DFP2 = (A[g][c] + 1 / a[g][c]) * (Q[:, c] / (mu[g][c] * (a[g][c] * mu[g][c] + 1)))
    DFP = DFP1 + DFP2

    DFQ1 = -P[g, :] / (a[g][c] * mu[g][c])
    DFQ2 = (A[g][c] + 1 / a[g][c]) * (P[g, :] / (mu[g][c] * (a[g][c] * mu[g][c] + 1)))
    DFQ = DFQ1 + DFQ2

    DFAg1 = (log(ag[g][0]) - 1) / (ag[g][0] ** 2)
    DFAg2 = log(mu[g][c]) / (ag[g][0] ** 2)
    DFAg3 = psi(1 / ag[g][0]) / (ag[g][0] ** 2)
    DFAg4 = -psi(A[g][c] + (1 / ag[g][0])) / (ag[g][0] ** 2)
    DFAg5 = log(1 + 1 / (ag[g][0] * mu[g][c])) / (ag[g][0] ** 2)
    DFAg6 = (A[g][c] + 1 / ag[g][0]) * (1 / (ag[g][0] * (mu[g][c] ** 2 + 1)))
    DFAg = DFAg1 + DFAg2 + DFAg3 + DFAg4 + DFAg5 + DFAg6

    DFAc1 = (log(ac[0][c]) - 1) / (ac[0][c] ** 2)
    DFAc2 = log(mu[g][c]) / (ac[0][c] ** 2)
    DFAc3 = psi(1 / ac[0][c]) / (ac[0][c] ** 2)
    DFAc4 = -psi(A[g][c] + (1 / ac[0][c])) / (ac[0][c] ** 2)
    DFAc5 = log(1 + 1 / (ac[0][c] * mu[g][c])) / (ac[0][c] ** 2)
    DFAc6 = (A[g][c] + 1 / ac[0][c]) * (1 / (ac[0][c] * (mu[g][c] ** 2 + 1)))
    DFAc = DFAc1 + DFAc2 + DFAc3 + DFAc4 + DFAc5 + DFAc6

    return DFP, DFQ, DFAg, DFAc