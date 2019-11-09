import torch
import pickle
import seaborn as sns
import numpy as np
import pandas as pd
from umap import UMAP
from sklearn.cluster import SpectralClustering
from scvi.inference import UnsupervisedTrainer
from scvi.models import VAE
save_path = '../CSF/Notebooks/'

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
%matplotlib inline

def ES_fast(score, s, p, interval):
    N = len(s)
    N_H = np.sum(s==1)
    m = 1/(N-N_H)
    power = np.abs(score)**p
    N_R = np.sum(power[s==1])
    h = power / N_R
    ES = [0]
    hit = 0
    miss = 0

    for i in np.arange(0, (len(power)-interval),interval):
        x = np.arange(i,i+interval,1)
        si = s[x]
        hit = hit + np.sum(h[x][si==1])
        miss = miss + m*np.sum(si==0)
        ES.append(hit-miss)

    return(ES)
