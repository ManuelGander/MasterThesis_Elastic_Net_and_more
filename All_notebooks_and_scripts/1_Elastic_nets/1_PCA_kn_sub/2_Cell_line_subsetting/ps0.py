import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py <argument>, wrong number of arguments")
    sys.exit(1)
i0 = int(sys.argv[1])
i1 = int(sys.argv[2])
i2 = int(sys.argv[3])
i3 = int(sys.argv[4])


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
import time
from os import listdir

import warnings
from scipy.stats import ConstantInputWarning
from sklearn.exceptions import ConvergenceWarning

# Ignore ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


alphas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
l1_ratios = [0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99]
seeds = [0,1,2,3]
frs = [0.25, 0.5, 0.75, 1.0]
kn=10**9
reps=100


alpha = alphas[i0]
l1_ratio = l1_ratios[i1]
seed = seeds[i2]
fr = frs[i3]

dataset = 'RNA'
source = 'GDSC2'


key = f'{alpha}_{l1_ratio}_{seed}_{fr}.pkl'
mypath='/home/icb/manuel.gander/Atl/data/dones_sub_new'
keys=listdir(mypath)
if not key in keys:


    D_prot, features, celllines = utils.load_dataset(dataset)
    dfv = utils.prep_viability_AUCs2(D_prot)
    dfv = dfv[dfv['Source']==source].copy()
    dfv = dfv[np.isfinite(dfv['AUC'])].copy()

    vc = dfv['PubChem_CID'].value_counts()
    drugs = vc[vc>150].index 
    dfv = dfv[dfv['PubChem_CID'].isin(drugs)].copy()
    ns=[]
    ds=[]
    for d,df in list(dfv[['PubChem_CID', 'Cello', 'sensitivity']].groupby('PubChem_CID')):
        ns.append((df['sensitivity']=='yes').sum())
        ds.append(d)
    dfd = pd.DataFrame({'PubChem_CID':ds, 'n_sensitive':ns})
    drugs = list(dfd[dfd['n_sensitive']>25]['PubChem_CID'])


    ccls = sorted(set(dfv['Cello']))
    np.random.seed(seed)
    ccls_sel = np.random.choice(ccls, int(len(ccls)*fr))
    dfv = dfv[dfv['Cello'].isin(ccls_sel)].copy()



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
    Df['source']=source
    Df['dataset']=dataset
    Df['alpha']=alpha
    Df['l1_ratio']=l1_ratio
    Df['seed']=seed
    Df['fr']=fr
    Df['kn']=kn

    Df.to_pickle(f'/home/icb/manuel.gander/Atl/data/dones_sub_new/{alpha}_{l1_ratio}_{seed}_{fr}.pkl')

