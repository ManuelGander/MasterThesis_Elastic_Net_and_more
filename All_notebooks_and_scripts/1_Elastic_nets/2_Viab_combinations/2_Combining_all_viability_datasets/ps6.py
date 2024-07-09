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

Path='/home/icb/manuel.gander/Atl/data'
def z_transform(M0):
    M0=M0-np.array(M0.mean(1))[:,None]
    M0=M0/np.array(M0.std(1))[:,None]
    return(M0)


alphas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
l1_ratios = [0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99]
datasets = ['atl_only_phos', 'atl_only_prot', 'atlantic', 'RNA', 'all', 'RNA_and_full', 'RNA_nz', 'RNA_nz_full']

kn = 10**6
reps=100

source = 'GDSC2'
alpha = alphas[i0]
l1_ratio = l1_ratios[i1]
dataset = datasets[i2]


def get_M(dataset):

    M0=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/lwn_imp_M0.pkl'))
    M1=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/als_imp_M1.pkl'))

    M0.columns=[a.split('_')[0] for a in M0.columns]
    M1.columns=[a.split('_')[0] for a in M1.columns]
    M0 = M0.loc[:,~M0.columns.duplicated()].copy()
    M1 = M1.loc[:,~M1.columns.duplicated()].copy()
    M2=z_transform(pd.read_pickle(f'{Path}/CCLE_RNA.pkl'))
    M2 = M2.loc[:,~M2.columns.duplicated()].copy()
    shared_cols = sorted(set(M0.columns)&set(M2.columns))
    M0=M0[shared_cols].copy()
    M1=M1[shared_cols].copy()
    M2=M2[shared_cols].copy()

    M3=pd.concat([M0,M1]).copy()
    M4=pd.concat([M3,M2]).copy()
    
    M5=pd.concat([M1,M2]).copy()
    
    
    M6=pd.read_pickle(f'{Path}/CCLE_RNA.pkl')
    M6 = M6.loc[:,~M6.columns.duplicated()].copy()
    M6=M6[shared_cols].copy()
    
    
    if dataset=='atl_only_phos':
        return(M0)
    elif dataset=='atl_only_prot':
        return(M1)
    elif dataset=='atlantic':
        return(M3)
    elif dataset=='RNA':
        return(M2)
    elif dataset=='all':
        return(M4)
    elif dataset=='RNA_and_full':
        return(M5)
    elif dataset=='RNA_nz':
        return(M6)
    else:
        print('Error: Unknown dataset')
        
        
M = get_M(dataset)
D_prot={}
for i,a in enumerate(M.columns):
    D_prot[a+'_i']=np.array(M[a])
features=list(M.index)
celllines=[a.split('_')[0] for a in D_prot.keys()]

with open('/home/icb/manuel.gander/Atl/data/common_drugs.pkl', 'rb') as file:
    Dd = pickle.load(file)
drugs = Dd[source]

dfv = utils.prep_viability_AUCs2(D_prot)
dfv = dfv[dfv['Source']==source].copy()

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


Df.to_pickle(f'/home/icb/manuel.gander/Atl/data/allf/{alpha}_{l1_ratio}_{dataset}.pkl')




