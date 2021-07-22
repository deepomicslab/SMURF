# _*_ coding: utf-8 _*_
"""
Time:     2021/7/22 17:33
Author:   WANG Bingchen
Version:  V 0.1
File:     lossFunc.py
Describe: 
"""
from scipy.special import gammaln
from numpy import log


def LossFunctionConstantVariance(mu, v, y):
    mu = max(mu, 0.01)
    v = max(v, 1e-09)
    fun1 = ((mu ** 2) / v) * log(mu)
    fun2 = -(mu ** 2 * (1 / v) * log(v))
    fun3 = -gammaln(mu ** 2 / v)
    fun4 = gammaln(y + mu**2 / v)
    fun5 = -(y + mu ** 2 / v) * log(1 + mu / v)
    fun = fun1 + fun2 + fun3 + fun4 + fun5
    return fun


def LossFunctionFano(mu, b, y):
    mu = max(mu, 0.01)
    b = max(b, 1e-09)
    fun1 = - (mu/b) * log(b)
    fun2 = -gammaln(mu/b)
    fun3 = gammaln(y + (mu/b))
    fun4 = -(y + (mu/b))*log(1 + (1/b))
    fun = fun1 + fun2 + fun3 + fun4
    return fun


def LossFunctionConstantCoefficientVariation(mu, a, y):
    mu = max(mu, 0.01)
    a = max(a, 1e-09)
    fun1 = -(1/a) * log(a)
    fun2 = -(1/a) * log(mu)
    fun3 = -gammaln(1/a)
    fun4 = gammaln(y + 1/a)
    fun5 = -(y + 1/a) * log(1 + 1/(a * mu))
    fun = fun1 + fun2 + fun3 + fun4 + fun5
    return fun
