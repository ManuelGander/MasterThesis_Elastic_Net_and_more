import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py <argument>, wrong number of arguments")
    sys.exit(1)
i0 = int(sys.argv[1])
i1 = int(sys.argv[2])
i2 = int(sys.argv[3])
i3 = int(sys.argv[4])


def z_transform(M0):
    M0=M0-np.array(M0.mean(1))[:,None]
    M0=M0/np.array(M0.std(1))[:,None]
    return(M0)

import pandas as pd
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


import warnings
from scipy.stats import ConstantInputWarning
from sklearn.exceptions import ConvergenceWarning

# Ignore ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


alphas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
l1_ratios = [0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99]
imputations = ['knn', 'als', 'hlm', 'lwn', 'zer', 'fals']
z_scores = [True, False]

kn=10**9
reps=100


alpha = alphas[i0]
l1_ratio = l1_ratios[i1]
data = 'Gygi'
imputation = imputations[i2]
z_score = z_scores[i3]
source = 'MR_NCI60'



df = pd.read_pickle(f'{imputation}_imp_M0.pkl')
    
if z_score:
    df = z_transform(df)

df.columns = [a+'_'+str(i) for i,a in enumerate(df.columns)]

D_prot={}
for a in df.columns:
    D_prot[a]=np.array(df[a])
features=list(df.index)
celllines=sorted(set([a.split('_')[0] for a in D_prot.keys()]))

dfv = utils.prep_viability_AUCs2(D_prot)
dfv = dfv[dfv['Source']==source].copy()
drugs = sorted(set(dfv['PubChem_CID']))



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
    if len(ccls)>30:
        D_holdout = utils.get_holdouts_balanced(dfvs, reps)
        if len(D_holdout)>0:
            for hs in range(reps):
                holdouts=D_holdout[hs]
                dfv0, dfh=utils.split_of_validation_cell_lines(dfvs, ccls=holdouts)
                train_input=utils.dfv_to_train_arrays(dfv0, D_prot)
                mean = train_input[1].mean()

                # Preselection of kn features and inference
                predictions, ys = utils.calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio)
                p=predictions[0]
                y=ys[0]

                res['pearsons'].append(scipy.stats.pearsonr(p, y).statistic)
                res['spearmans'].append(scipy.stats.spearmanr(p, y).statistic)
                res['RMSE'].append(np.sqrt(np.mean((p-y)**2)))
                res['l1_ratios'].append(np.mean(abs(p-y))/np.mean(abs(mean-y)))
                res['RMSE_mean_model'].append(np.sqrt(np.mean((mean-y)**2)))

    df = pd.DataFrame(res)
    df['PubChem_CID']=dr
    df['ind']=range(len(df))
    df['n_ccls'] = len(dfvs)
    frames.append(df)
Df = pd.concat(frames, ignore_index=True)

Df['imputation']=imputation
Df['z_score']=z_score
Df['source']=source
Df['dataset']=data
Df['alpha']=alpha
Df['l1_ratio']=l1_ratio


Df.to_pickle(f'/home/icb/manuel.gander/Atl/data/dones_imp_Gygi/{alpha}_{l1_ratio}_{imputation}_{data}_{z_score}.pkl')

