import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def ES_fast(score, s, p, interval):
    N = len(s)
    N_H = np.sum(s == 1)
    m = 1 / (N - N_H)
    power = np.abs(score) ** p
    N_R = np.sum(power[s == 1])
    h = power / N_R
    ES = [0]
    hit = 0
    miss = 0

    for i in np.arange(0, (len(power) - interval), interval):
        x = np.arange(i, i + interval, 1)
        si = s[x]
        hit = hit + np.sum(h[x][si == 1])
        miss = miss + m * np.sum(si == 0)
        ES.append(hit - miss)

    return ES


def matchedlist(norm_X, genelist, genenames, ngenes):
    mean_exprs = np.asarray(norm_X[:, :].mean(axis=0)).ravel()
    plt.hist(np.log10(1 + mean_exprs))
    matched = []
    for x in genelist:
        if x in genenames:
            exprs = np.mean(norm_X[:, genenames == x])
            plt.axvline(np.log10(1 + exprs), color="r")
            diff = np.abs(mean_exprs - exprs)
            idx = np.argsort(np.argsort(diff))
            idx = np.where(idx <= ngenes)[0]
            res = genenames[idx]
            res = res[res != x]
            matched.append(res)
            exprs = np.mean(norm_X[:, idx])
            plt.axvline(np.log10(1 + exprs), color="b")
    return matched
