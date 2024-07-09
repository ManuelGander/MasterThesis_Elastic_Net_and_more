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
import scipy
import time
import numpy as np
import pickle

# torch
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import ExponentialLR
from tqdm import tqdm
from sklearn.decomposition import PCA
import sys
sys.path.append('/home/icb/manuel.gander/Atl/notebooks/')
import utils
import importlib
utils = importlib.reload(utils)
device = 'cuda'
Path = '/home/icb/manuel.gander/Atl/data'



# 3 3 2 3
datasets = ['RNA', 'identity', 'atl_only_full', 'atl_only_phos']
dataset = datasets[i0]

n_componentss = [2,5,10,30]
n_components = n_componentss[i1]

methods = ['als', 'kino', 'rdkit']
method = methods[i2]

dense_dims = [5,10,20,50]
dense_dim = dense_dims[i3]


dfv_train = pd.read_pickle(f'{dataset}_dfv_train.pkl')
dfv_test = pd.read_pickle(f'{dataset}_dfv_test.pkl')

Za = pd.read_pickle(f'{dataset}_ALS_drug_embedding.pkl')
Dai = dict(zip([str(a) for a in Za.index], range(len(Za))))
Zan = Za.values
Zk = pd.read_pickle(f'{dataset}_Kinobead_drug_embedding.pkl')
Dki = dict(zip([str(a) for a in Zk.index], range(len(Zk))))
Zkn = Zk.values
Zr = pd.read_pickle(f'{dataset}_rdkit_drug_embedding.pkl')
Dri = dict(zip([str(a) for a in Zr.index], range(len(Zr))))
Zrn = Zr.values

# only keep the drugs that are shared between all of them to compare it all properly
drugs = sorted([str(a) for a in set(Dai.keys())&set(Dki.keys())&set(Dri.keys())])

dfv_train = dfv_train[dfv_train['PubChem_CID'].isin(drugs)].copy()
dfv_test = dfv_test[dfv_test['PubChem_CID'].isin(drugs)].copy()
dfv_train_mean = dfv_train[['PubChem_CID', 'AUC']].groupby(['PubChem_CID']).mean().reset_index().copy()
Dvm = dict(zip(dfv_train_mean.PubChem_CID, dfv_train_mean.AUC))
dfv_test['mean_AUC'] = dfv_test['PubChem_CID'].map(Dvm)

# get torch arrays of feature data
with open(f'/home/icb/manuel.gander/Atl/data/pcas/{dataset}_{n_components}.pkl', 'rb') as file:
    D_prot2, features, cellos = pickle.load(file)
P = pd.DataFrame(np.vstack(list(D_prot2.values())), index=[a.split('_')[0] for a in D_prot2.keys()])
Dpi = dict(zip(P.index, range(len(P))))
Pn = P.values


if method == 'als':
    Di = Dai
    Z = Zan
elif method == 'kino':
    Di = Dki
    Z = Zkn
elif method == 'rdkit':
    Di = Dri
    Z = Zrn
    
D = Z[dfv_train['PubChem_CID'].map(Di),:].copy()
C = Pn[dfv_train['Cello'].map(Dpi),:]

D2 = Z[dfv_test['PubChem_CID'].map(Di),:].copy()
C2 = Pn[dfv_test['Cello'].map(Dpi),:]


Dt = torch.tensor(D, device=device, dtype=torch.float32)
Ct = torch.tensor(C, device=device, dtype=torch.float32)
vt = torch.tensor(dfv_train['AUC'].values, device=device, dtype=torch.float32)

D2t = torch.tensor(D2, device=device, dtype=torch.float32)
C2t = torch.tensor(C2, device=device, dtype=torch.float32)
v2t = torch.tensor(dfv_test['AUC'].values, device=device, dtype=torch.float32)
v2mt = torch.tensor(dfv_test['mean_AUC'].values, device=device, dtype=torch.float32)


import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.nn.init as init

class Net(nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        self.d1 = nn.Linear(Ct.shape[1], dense_dim)
        self.d2 = nn.Linear(Dt.shape[1], dense_dim)
        self.d3 = nn.Linear(dense_dim*2, dense_dim)
        self.d4 = nn.Linear(dense_dim, dense_dim)
        self.d5 = nn.Linear(dense_dim, 1)
        
    def forward(self, x, y):
        x = F.relu(self.d1(x))
        y = F.relu(self.d2(y))
        z = torch.concat([x,y], axis=1)
        z = F.relu(self.d3(z))
        z = F.relu(self.d4(z))
        z = self.d5(z)
        #z = 1-F.sigmoid(self.d5(z))
        return(z)
net = Net().to(device)
print(net)
params = list(net.parameters())


mean_model_loss = torch.norm(v2t-v2mt)

tr_losses = []
test_losses = []

lr0=1e-2
lr1=1e-6
steps=5*10**5
gamma=np.exp(np.log(lr1/lr0)/steps)


for i in range(steps):
    lr=lr0*gamma**i
    
    out = net(Ct, Dt).flatten()
    loss = torch.norm(out-vt)
    your_optimizer = optim.Adam(list(net.parameters()), lr=lr)
    your_optimizer.zero_grad()
    loss.backward()
    your_optimizer.step()
    
    
    if i%10**3==0:
        tr_losses.append(float(loss))
        out2 = net(C2t, D2t).flatten()
        test_loss = torch.norm(out2-v2t)
        print(float(test_loss/mean_model_loss))
        test_losses.append(float(test_loss))
        

dfs = pd.DataFrame({'dataset':dataset, 'n_components':n_components, 'method':method, 'dense_dim':dense_dim,
             'tr_losses':tr_losses, 'test_losses':test_losses, 'mean_loss':float(mean_model_loss)})
dfs.to_csv(f'{Path}/dr_NN/{dataset}_{n_components}_{method}_{dense_dim}.csv')
