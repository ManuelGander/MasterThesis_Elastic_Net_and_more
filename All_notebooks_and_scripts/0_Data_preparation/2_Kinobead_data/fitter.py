import pandas as pd
import numpy as np
import csv
import pickle
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit
from tqdm import tqdm
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import curve_curator
from curve_curator.models import _Model, MeanModel, LogisticModel
from curve_curator.quantification import fit_model

M0 = MeanModel()
M1 = LogisticModel()
fit_params = {'type': 'OLS', 'speed':'standard', 'control_fold_change':False}
# I took these as standart parameters https://github.com/kusterlab/curve_curator/blob/main/example_toml_files/all_parameters.toml
f_statistic_params = {'optimized_dofs':False, 'scale':1, 'loc':0.12}

def process_dfs(dfs):
    key = dfs['key'].iloc[0]
    xs = list(dfs['log_dose'])
    ys = list(dfs['Viability'])

    xs.append(10)
    ys.append(0)
    
    ns=round(dfs['Dose'].value_counts().median())
    xs = np.array([*xs, *[min(xs)-3]*ns])
    ys = np.array([*ys, *[1]*ns])

    res = fit_model(pd.Series(ys), xs, M0, M1, fit_params, f_statistic_params)
    #params = res[:4]
    #rmse, r2 = res[5:7]
    return(key, res)
