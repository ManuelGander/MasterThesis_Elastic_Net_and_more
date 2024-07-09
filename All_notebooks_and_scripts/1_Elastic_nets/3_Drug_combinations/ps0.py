import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 5:
    print("Usage: python script.py <argument>, wrong number of arguments")
    sys.exit(1)
i0 = int(sys.argv[1])
i1 = int(sys.argv[2])
i2 = int(sys.argv[3])
i3 = int(sys.argv[4])

# Check if there are keys missing
def remove_non_alphanumeric(input_string):
    # Muster für alle alphanumerischen Zeichen und Leerzeichen
    pattern = re.compile(r'[^a-zA-Z0-9 ]+')
    # Ersetze alle nicht-alphanumerischen Zeichen mit einem leeren String
    result = pattern.sub('', input_string)
    return result


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
import pickle
import re

import warnings
from scipy.stats import ConstantInputWarning
from sklearn.exceptions import ConvergenceWarning

# Ignore ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)

datasets = ['atlantic', 'atl_only_phos', 'atl_only_full', 'Gygi', 'kinase_scores', 'RNA', 'drug_scores', 'identity']
alphas = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
l1_ratios = [0.01, 0.03, 0.1, 0.3, 0.5, 0.9, 0.99]
classes = ['Apoptosis signaling',
 'Cdk inhibitor',
 'Checkpoint signaling',
 'Cytoskeleton',
 'DNA degradation',
 'DNA repair',
 'Heat / shock',
 'Histone acetyltransferase',
 'Histone methyltransferase',
 'Inflammation',
 'Inositole phosphorylation',
 'MAPK signaling',
 'Nfkg signaling',
 'RTK',
 'Topoisomerase',
 'Transcription',
 'mTOR signaling']

kn=10**9
reps=100

alpha = alphas[i0]
l1_ratio = l1_ratios[i1]
dataset = datasets[i2]
c = classes[i3]
cs = remove_non_alphanumeric(c)
source = 'GDSC2'


Path = '/home/icb/manuel.gander/Atl/data'
dfd = pd.read_csv(f'{Path}/drugs_annotated.txt', sep='\t')
dfd['Drug'] = [a[:-1] for a in dfd['Drug']]
dfd.index=range(len(dfd))
dfd['dr'] = [remove_non_alphanumeric(a).lower() for a in dfd['Drug']]
D_cl_to_drugs = {c:list(dfd[dfd['MOA']==c]['dr']) for c in classes}

D_prot, features, celllines = utils.load_dataset(dataset)
dfv = utils.prep_viability_AUCs2(D_prot)
dfv = dfv[dfv['Source']==source].copy()
dfv['dr'] = [remove_non_alphanumeric(a).lower() for a in dfv['Dr_repr_name']]

with open('/home/icb/manuel.gander/Atl/data/common_drugs.pkl', 'rb') as file:
    Dd = pickle.load(file)
drugs = Dd[source]






frames=[]
dfvs = dfv[dfv['dr'].isin(D_cl_to_drugs[c])].copy()
dfvs['PubChem_CID_initially'] = dfvs['PubChem_CID'].copy()
dfvs['PubChem_CID'] = c

ccls=sorted(set(dfvs['Cello']))
if len(ccls)>30:
    D_holdout = utils.get_holdouts_balanced(dfvs, reps, source0='all')
    if len(D_holdout)>0:
        for hs in tqdm(range(reps)):
            holdouts=D_holdout[hs]
            dfv0, dfh=utils.split_of_validation_cell_lines(dfvs, ccls=holdouts)
            train_input=utils.dfv_to_train_arrays(dfv0, D_prot)
            mean = train_input[1].mean()

            # Preselection of kn features and inference
            predictions, ys = utils.calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio)
            p=predictions[0]
            y=ys[0]
            dfh['pred'] = p
            drugs = sorted(set(dfh['PubChem_CID_initially']))
            for d in drugs:
                dfhs = dfh[dfh['PubChem_CID_initially']==d].copy()
                if len(dfhs)>4:
                    mean = dfv0[dfv0['PubChem_CID_initially']==d]['AUC'].mean()


                    pc = scipy.stats.pearsonr(dfhs['AUC'], dfhs['pred']).statistic
                    rmse = np.sqrt(np.mean((dfhs['AUC']-dfhs['pred'])**2))
                    rmse_mean = np.sqrt(np.mean((dfhs['AUC']-mean)**2))

                    df = pd.DataFrame({'PubChem_CID':[d], 'Class':c, 'pearsons':pc, 'RMSE':rmse, 'RMSE_mean_model':rmse_mean,
                                  'ind':hs, 'n_ccls':len(dfvs[dfvs['PubChem_CID_initially']==d])})
                    frames.append(df)
        Df = pd.concat(frames, ignore_index=True)
        Df['alpha'] = alpha
        Df['l1_ratio'] = l1_ratio
        Df['dataset'] = dataset
        Df['source']=source
        Df['kn'] = kn
            
        
        Df.to_pickle(f'/home/icb/manuel.gander/Atl/data/dr_cl_comb2/{alpha}_{l1_ratio}_{dataset}_{cs}.pkl')
    else:
        pd.DataFrame([]).to_pickle(f'/home/icb/manuel.gander/Atl/data/dr_cl_comb2/{alpha}_{l1_ratio}_{dataset}_{cs}.pkl')
else:
    pd.DataFrame([]).to_pickle(f'/home/icb/manuel.gander/Atl/data/dr_cl_comb2/{alpha}_{l1_ratio}_{dataset}_{cs}.pkl')
