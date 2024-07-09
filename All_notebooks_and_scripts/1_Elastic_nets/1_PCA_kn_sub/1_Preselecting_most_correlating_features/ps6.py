import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py <argument>, wrong number of arguments")
    sys.exit(1)
i0 = int(sys.argv[1])
i1 = int(sys.argv[2])
i2 = int(sys.argv[3])


import pandas as pd
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import importlib
from tqdm import tqdm
import sys
sys.path.append('/home/icb/manuel.gander/Atl/notebooks/')
import utils
utils = importlib.reload(utils)
import scipy
import pickle

import warnings
from scipy.stats import ConstantInputWarning
from sklearn.exceptions import ConvergenceWarning

# Ignore ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


alphas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
l1_ratios = [0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99]
kns = [10, 30, 100, 300, 1000, 3000, 10000, 20000]
reps=100

alpha = alphas[i0]
l1_ratio = l1_ratios[i1]
kn = kns[i2]

source = 'MR_NCI60'
#source = 'GDSC2'
dataset = 'atl_only_phos'


D_prot, features, celllines = utils.load_dataset(dataset)
dfv = utils.prep_viability_AUCs2(D_prot)
dfv = dfv[dfv['Source']==source].copy()
dfv = dfv[np.isfinite(dfv['AUC'])].copy()    

with open('/home/icb/manuel.gander/Atl/data/common_drugs.pkl', 'rb') as file:
    Dd = pickle.load(file)
drugs = Dd[source]


frames=[]
for dr in tqdm(drugs):
    dfvs=dfv[dfv['PubChem_CID']==dr].copy()

    ccls=sorted(set(dfvs['Cello']))
    res={}
    res['pearsons']=[]
    res['spearmans']=[]
    res['l1_ratios']=[]
    res['RMSE']=[]
    res['RMSE_mean_model']=[]
    res['td']=[]
    if len(ccls)>30:
        D_holdout = utils.get_holdouts_balanced(dfvs, reps, source0='all')
        if len(D_holdout)>0:
            for hs in range(reps):
                holdouts=D_holdout[hs]
                dfv0, dfh=utils.split_of_validation_cell_lines(dfvs, ccls=holdouts)
                train_input=utils.dfv_to_train_arrays(dfv0, D_prot)
                mean = train_input[1].mean()

                # Preselection of kn features and inference
                t0=time.time()
                predictions, ys = utils.calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio)
                td = time.time()-t0
                p=predictions[0]
                y=ys[0]

                res['pearsons'].append(scipy.stats.pearsonr(p, y).statistic)
                res['spearmans'].append(scipy.stats.spearmanr(p, y).statistic)
                res['RMSE'].append(np.sqrt(np.mean((p-y)**2)))
                res['l1_ratios'].append(np.mean(abs(p-y))/np.mean(abs(mean-y)))
                res['RMSE_mean_model'].append(np.sqrt(np.mean((mean-y)**2)))
                res['td'].append(td)

    df = pd.DataFrame(res)
    df['PubChem_CID']=dr
    df['ind']=range(len(df))
    df['n_ccls'] = len(dfvs)
    frames.append(df)
Df = pd.concat(frames, ignore_index=True)
Df['source']=source
Df['dataset']=dataset
Df['alpha']=alpha
Df['l1_ratio']=l1_ratio
Df['kn'] = kn


Df.to_pickle(f'/home/icb/manuel.gander/Atl/data/dones_kn2/{alpha}_{l1_ratio}_{kn}.pkl')




