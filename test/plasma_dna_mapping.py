import numpy as np
from cvxopt import solvers, matrix
from tqdm import tqdm
from scipy import optimize


def QP(reference, mixtureData):
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    s = np.rank(reference)
    for i in range(numOfSamples):
        P = matrix(np.dot(np.transpose(reference), reference))
        Q = matrix(-1 * np.dot(np.transpose(reference), mixtureData[:, i]))
        G = matrix(-1 * np.eye(tissueNum))
        H = matrix(np.zeros([tissueNum, 1]))
        A = matrix(np.ones([1, tissueNum]))
        B = matrix(np.ones([1, 1]))
        solvers.options["show_progress"] = False
        sol = solvers.qp(P, Q, G, H, A, B)
        t = np.asarray(sol["x"])
        proportionDeconvolution[:, i] = t[:, 0]
    return proportionDeconvolution


def PlasmaDNATissueMapping(reference, mixture):
    newMarkerNumber = np.size(reference, 0)
    cellTypeNumber = np.size(reference, 1)
    selectedCpG = np.arange(newMarkerNumber)
    area = np.zeros(newMarkerNumber)
    area = area.astype(np.bool)
    for i in tqdm(range(newMarkerNumber)):
        temp = reference[i, :]
        SD = np.std(temp)
        methylationMax = np.max(temp) / np.min(temp)
        variability = SD / np.mean(temp)
        if variability >= 0.25 and methylationMax > 0.2:
            area[i] = True
        for j in range(cellTypeNumber):
            if temp[j] < np.mean(temp) - 3 * SD or temp[j] > np.mean(temp) + 3 * SD:
                area[i] = True
    selectedCpG = selectedCpG[area]
    referenceNew = reference[selectedCpG, :]
    mixtureNew = mixture[selectedCpG, :]
    numberOfSamples = np.size(mixtureNew, 1)
    proportionDeconvolution = QP(referenceNew, mixtureNew)
    return proportionDeconvolution


def NNLS(reference, mixtureData):
    tissueNum = np.size(reference, 1)
    numOfSamples = np.size(mixtureData, 1)
    proportionDeconvolution = np.zeros([tissueNum, numOfSamples])
    for i in tqdm(range(numOfSamples)):
        mixture = optimize.nnls(reference, mixtureData[:, i])
        test = mixture[0]
        t = test / np.sum(test)
        proportionDeconvolution[:, i] = t
    return proportionDeconvolution
