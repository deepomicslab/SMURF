# _*_ coding: utf-8 _*_
"""
Time:     2021/7/22 17:29
Author:   WANG Bingchen
Version:  V 0.1
File:     scend.py
Describe: 
"""

from .lossFunc import *
from .initialParams import *
import pandas as pd
from .optimize import CVOptimize, FanoOptimize, CCVOptimize
from . import utiles
import warnings
warnings.filterwarnings("ignore")


class SCEnd():
    def __init__(self, n_features=20, steps=10, alpha=1e-5, eps=10, noise_model="Fano", normalize=True, calculateIntialNoiseFactor=False, estimate_only=False):

        self.K = n_features
        self.batchSize = n_features * 10
        self.steps = steps
        self.alpha = alpha
        self.eps = eps
        self.normalize = normalize


        self.noise_model = noise_model
        self.calculateLossFunc = True
        self.estmate_only = estimate_only
        self.calculateIntialNoiseFactor = calculateIntialNoiseFactor

    def _check_params(self):
        """Check parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as n='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utiles.check_positive(n_features=self.K, eps=self.eps, batchsize=self.batchSize, alpha=self.alpha,
                              steps=self.steps)
        utiles.check_int(n_features=self.K, steps=self.steps)
        utiles.check_between(v_min=0, v_max=min(self.genes, self.cells), n_features=self.K)
        utiles.check_bool(normalize=self.normalize)
        utiles.check_bool(iteration=self.calculateIntialNoiseFactor)
        utiles.check_noise_model(noise_model=self.noise_model)

    def MFConstantVariance(self):
        P = self.P
        Q = self.Q
        A = self.A
        mu = np.dot(P, Q)
        Vc = np.ones((1, self.cells)) * 1.0

        if self.calculateIntialNoiseFactor:
            Vg = getv(self.A, mu)
            v = np.dot(Vg, Vc)
        else:
            Vg = np.ones((self.genes, 1))
            v = np.dot(Vg, Vc)

        LossFuncPre = 0
        LossFunc = 0
        Lossf = []

        nonZerosInd = np.nonzero(A)
        numNonZeros = nonZerosInd[0].size
        nonZeroElem = []

        for indice in range(numNonZeros):
            g = nonZerosInd[0][indice]
            c = nonZerosInd[1][indice]
            element = [g, c]
            nonZeroElem.append(element)
        nonZeroElem = np.array(nonZeroElem)
        print("start iteration")


        for step in range(self.steps):

            if self.calculateLossFunc:
                for element in nonZeroElem:
                    g = element[0]
                    c = element[1]

                    LossFuncCG = LossFunctionConstantVariance(mu[g][c], v[g][c], A[g][c])
                    LossFunc = LossFunc + LossFuncCG
            else:
                LossFunc = (step+1)*(self.eps+1)


            if abs(LossFunc - LossFuncPre) < self.eps:
                Lossf.append(LossFunc)
                print("already converge")
                break
            else:
                Lossf.append(LossFunc)
                LossFuncPre = LossFunc
                LossFunc = 0

                np.random.shuffle(nonZeroElem)

                batchElements = nonZeroElem[0:self.batchSize - 1, :]

                for element in batchElements:
                    g = element[0]
                    c = element[1]
                    if mu[g][c] <= 0.01:
                        mu[g][c] = 0.01
                    if Vg[g][0] <= 1e-09:
                        Vg[g][0] = 1e-09

                    dP, dQ, dVg = CVOptimize(A, P, Q, mu, Vg, Vc, v, g, c)
                    P[g, :] = P[g, :] + self.alpha * dP
                    Q[:, c] = Q[:, c] + self.alpha * dQ
                    Vg[g, :] = Vg[g, :] + self.alpha * dVg

                mu = np.dot(P, Q)
                v = np.dot(Vg, Vc)

            print("number of iteration: ", step+1,"/",self.steps)
        mu = np.dot(P, Q)
        mu[mu < 0] = 0
        return mu, P, Q


    def MFFano(self):
        P = self.P
        Q = self.Q
        A = self.A
        mu = np.dot(P, Q)
        bc = np.ones((1, self.cells)) * 1.0

        if self.calculateIntialNoiseFactor:
            bg = getb(self.A, mu)
        else:
            bg = np.ones((self.genes, 1)) * 1.0

        b = np.dot(bg, bc)

        LossFunc = 0
        LossFuncPre = 0
        Lossf = []


        nonZerosInd = np.nonzero(A)
        numNonZeros = nonZerosInd[0].size
        nonZeroElem = []

        for indice in range(numNonZeros):
            g = nonZerosInd[0][indice]
            c = nonZerosInd[1][indice]
            element = [g, c]
            nonZeroElem.append(element)
        nonZeroElem = np.array(nonZeroElem)

        for step in range(self.steps):
            if self.calculateLossFunc:
                for element in nonZeroElem:
                    g = element[0]
                    c = element[1]

                    LossFuncCG = LossFunctionFano(mu[g][c], b[g][c], A[g][c])
                    LossFunc = LossFunc + LossFuncCG
            else:
                LossFunc = (step + 1)*self.eps


            if abs(LossFunc - LossFuncPre) < self.eps:
                Lossf.append(LossFunc)
                print("already converge")
                break
            else:
                Lossf.append(LossFunc)
                LossFuncPre = LossFunc
                LossFunc = 0

                np.random.shuffle(nonZeroElem)

                batchElements = nonZeroElem[0:self.batchSize - 1, :]
                for element in batchElements:
                    g = element[0]
                    c = element[1]
                    if mu[g][c] <= 0.01:
                        mu[g][c] = 0.01
                    if bg[g][0] <= 1e-09:
                        bg[g][0] = 1e-09
                    b[g][c] = np.dot(bg[g, :], bc[:, c])

                    dP, dQ, dbg = FanoOptimize(A, P, Q, mu, bg, bc, b, g, c)
                    P[g, :] = P[g, :] + self.alpha * dP
                    Q[:, c] = Q[:, c] + self.alpha * dQ
                    bg[g, :] = bg[g, :] + self.alpha * dbg
                mu = np.dot(P, Q)
                b = np.dot(bg, bc)
            self.alpha = self.alpha*(1 - (np.float(step)/np.float(self.steps)))
            print("number of iteration: ", step+1, "/", self.steps)
        mu = np.dot(P, Q)
        mu[mu < 0] = 0

        return mu, P, Q

    def MFConstCoeffiVariation(self):
        P = self.P
        Q = self.Q
        A = self.A
        mu = np.dot(P, Q)
        ac = np.ones((1, self.cells))
        if self.calculateIntialNoiseFactor:
            ag = geta(self.A, mu)
        else:
            ag = np.ones((self.genes, 1))

        a = np.dot(ag, ac)

        LossFunc = 0
        LossFuncPre = 0
        Lossf = []

        nonZerosInd = np.nonzero(A)
        numNonZeros = nonZerosInd[0].size
        nonZeroElem = []

        for indice in range(numNonZeros):
            g = nonZerosInd[0][indice]
            c = nonZerosInd[1][indice]
            element = [g, c]
            nonZeroElem.append(element)
        nonZeroElem = np.array(nonZeroElem)

        for step in range(self.steps):
            if self.calculateLossFunc:
                for element in nonZeroElem:
                    g = element[0]
                    c = element[1]
                    LossFuncCG = LossFunctionFano(mu[g][c], a[g][c], A[g][c])
                    LossFunc = LossFunc + LossFuncCG
            else:
                LossFunc = (step + 1)*self.eps


            if abs(LossFunc - LossFuncPre) < self.eps:
                Lossf.append(LossFunc)
                print("already converge")
                break
            else:
                Lossf.append(LossFunc)
                LossFuncPre = LossFunc
                LossFunc = 0

                np.random.shuffle(nonZeroElem)

                batchElements = nonZeroElem[0:self.batchSize - 1, :]
                for element in batchElements:
                    g = element[0]
                    c = element[1]
                    if mu[g][c] <= 0.01:
                        mu[g][c] = 0.01
                    if ag[g][0] <= 1e-09:
                        ag[g][0] = 1e-09
                    if ac[0][c] <= 1e-09:
                        ac[0][c] = 1e-09

                    dP, dQ, dag, dac = CCVOptimize(A, P, Q, mu, ag, ac, a, g, c)

                    P[g, :] = P[g, :] + self.alpha * dP
                    Q[:, c] = Q[:, c] + self.alpha * dQ
                    ag[g, :] = ag[g, :] + self.alpha * dag
                    ac[:, c] = ac[:, c] + self.alpha * dac
                mu = np.dot(P, Q)
                a = np.dot(ag, ac)
            print("number of iteration: ", step + 1, "/", self.steps)
        mu = np.dot(P, Q)
        mu[mu < 0] = 0

        return mu, P, Q

    def scend_impute(self, initialDataFrame):
        self.initialDataFrame = initialDataFrame
        self.genes = initialDataFrame.shape[0]
        self.cells = initialDataFrame.shape[1]
        self.genesNames = initialDataFrame._stat_axis.values.tolist()
        self.cellsNames = initialDataFrame.columns.values.tolist()


        print("Running SCEnd on {} cells and {} genes".format(self.cells, self.genes))



        if self.normalize:
            print("normalizing data by library size...")
            normalizedDataframe, self.size_factors = utiles.dataNormalization(self.initialDataFrame)
            self.A = normalizedDataframe.values
            self.P, self.Q = initialMatrices(normalizedDataframe.values, self.K)

            self._check_params()

            print("preprocessing data...")

            if self.noise_model == "CV":
                mu, P, Q = self.MFConstantVariance()

            if self.noise_model == "Fano":
                mu, P, Q = self.MFFano()

            if self.noise_model == "CCV":
                mu, P, Q = self.MFConstCoeffiVariation()




            newDataFrame = pd.DataFrame(mu, index=self.genesNames, columns=self.cellsNames)
            newDataFrame = newDataFrame * self.size_factors


            res = {}

            res["estimate"] = newDataFrame
            res["gene latent factor matrix"] = pd.DataFrame(P, index=self.genesNames, columns=None)
            res["cell latent factor matrix"] = pd.DataFrame(Q, index=None, columns=self.cellsNames)

            if self.estmate_only:
                return res["estimate"]
            else:
                return res

        else:
            self.A = self.initialDataFrame.values
            self.P, self.Q = initialMatrices(self.initialDataFrame.values, self.K)

            self._check_params()

            print("preprocessing data...")

            if self.noise_model == "CV":
                mu, P, Q = self.MFConstantVariance()

            if self.noise_model == "Fano":
                mu, P, Q = self.MFFano()

            if self.noise_model == "CCV":
                mu, P, Q = self.MFConstCoeffiVariation()

            newDataFrame = pd.DataFrame(mu, index=self.genesNames, columns=self.cellsNames)

            res = {}

            res["estimate"] = newDataFrame
            res["gene latent factor matrix"] = pd.DataFrame(P, index=self.genesNames, columns=None)
            res["cell latent factor matrix"] = pd.DataFrame(Q, index=None, columns=self.cellsNames)

            if self.estmate_only:
                return res["estimate"]
            else:
                return res
















