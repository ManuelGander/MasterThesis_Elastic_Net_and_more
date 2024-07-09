import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.linear_model import ElasticNet
import scipy

Path='/home/icb/manuel.gander/Atl/data'
def paerson_nb(x, y):
    A=x
    B=y
    A=A-A.mean(0)
    B=B-B.mean(0)
    B=np.divide(B,np.std(B,0)+10**-10)
    A=np.divide(A,np.std(A,0)+10**-10)
    A=np.nan_to_num(A,0)
    np.matmul(A.T,B)/A.shape[0]
    return np.matmul(A.T,B)/A.shape[0]

def split_of_validation_cell_lines(dfv, ccls):
    dfv0=dfv[~dfv['Cello'].isin(ccls)].copy()
    dfv1=dfv[dfv['Cello'].isin(ccls)].copy()
    return(dfv0, dfv1)

def dfv_to_train_arrays(dft, D_prot):
    X=np.vstack(dft['Cello_i'].map(D_prot))
    y=np.array(dft['AUC'])
    return(X,y)

def get_rs(iis, train_input):
    X=train_input[0][:,iis]
    y=train_input[1]
    model = LinearRegression()
    model.fit(X, y)
    predictions, y, mean = get_predictions(model, train_input, iis)
    rs=1-np.sum((y-predictions)**2)/np.sum((y-np.mean(y))**2)
    return float(rs)

def get_predictions(model, train_input, iis):
    X=train_input[0][:,iis]
    predictions = model.predict(X)
    y = train_input[1]
    mean = np.mean(y)
    return predictions, y, mean


def min_imputation(M0):
    mins=M0.min(axis=1)
    for i,j in zip(np.where(np.isnan(M0))[0], np.where(np.isnan(M0))[1]):
        M0.iloc[i,j]=mins[i]
    return (M0)

def z_transform(M0):
    M0=M0-np.array(M0.mean(1))[:,None]
    M0=M0/np.array(M0.std(1))[:,None]
    return(M0)


# This is the function I got from my friend that he uses for impuation of missing values in proteomics (LC-MS)

def impute_normal_down_shift_distribution(unimputerd_dataframe:pd.DataFrame ,column_wise=True, width=0.3, downshift=1.8, seed=2):
    """ 
    Performs imputation across a matrix columnswise
    https://rdrr.io/github/jdreyf/jdcbioinfo/man/impute_normal.html#google_vignette
    :width: Scale factor for the standard deviation of imputed distribution relative to the sample standard deviation.
    :downshift: Down-shifted the mean of imputed distribution from the sample mean, in units of sample standard deviation.
    :seed: Random seed
    
    """
    unimputerd_df = unimputerd_dataframe

    unimputerd_matrix = unimputerd_df.replace({pd.NA: np.nan}, inplace=True)
    
    unimputerd_matrix = unimputerd_df.to_numpy()
    columns_names = unimputerd_df.columns
    rownames = unimputerd_df.index
    unimputerd_matrix[~np.isfinite(unimputerd_matrix)] = None
    main_mean = np.nanmean(unimputerd_matrix)
    main_std = np.nanstd(unimputerd_matrix)
    np.random.seed(seed = seed)
    def impute_normal_per_vector(temp:np.ndarray,width=width, downshift=downshift):
        """ Performs imputation for a single vector """
        if column_wise:
            temp_sd = np.nanstd(temp)
            temp_mean = np.nanmean(temp)
        else:
            # over all matrix
            temp_sd = main_std
            temp_mean = main_mean

        shrinked_sd = width * temp_sd
        downshifted_mean = temp_mean - (downshift * temp_sd) 
        n_missing = np.count_nonzero(np.isnan(temp))
        temp[np.isnan(temp)] = np.random.normal(loc=downshifted_mean, scale=shrinked_sd, size=n_missing)
        if n_missing > 0:
            print 
        return temp
    final_matrix = np.apply_along_axis(impute_normal_per_vector, 0, unimputerd_matrix)
    final_df = pd.DataFrame(final_matrix)
    final_df.index = rownames
    final_df.columns = columns_names

    final_df = pd.concat([unimputerd_dataframe,final_df], axis=1) 
    
    return final_df

def purge_and_capitalize(ccl):
    ccl_str=ccl.replace('-', '')
    ccl_str=ccl_str.replace(' ', '')
    return(ccl_str.upper())

def map_celllines_to_cellosaurus(ccls):
    # Let's remove '-' and spaces and capitalize everything
    D_ccl={}

    for ccl in ccls:
        D_ccl[ccl]=purge_and_capitalize(ccl)
    ccl_red=sorted(set([D_ccl[a] for a in ccls]))
    len(ccl_red)
    

    dfc=pd.read_pickle(f'{Path}/Cellosaurus_CCLs.pkl')
    dfc=dfc[dfc['Species']=='NCBI_TaxID=9606; ! Homo sapiens (Human)']

    # Use the synonyms to map to the Cellosaurus ID
    cellosaurus_ids=list(dfc['Cellosaurus_ID'])
    syns=list(dfc['Synonyms'])

    syn_l=[]
    for s in syns:
        if '; ' in s:
            syn_l.append(s.split('; '))
        else:
            syn_l.append(s)
    D_cello={cellosaurus_ids[i]:cellosaurus_ids[i] for i in range(len(cellosaurus_ids))}

    for i in range(len(syn_l)):
        for s in syn_l[i]:
            D_cello[s]=cellosaurus_ids[i]
            
    ccl_red=sorted(set(ccl_red))
    for ccl in ccl_red:
        if ccl in cellosaurus_ids:
            D_cello[ccl]=ccl
            
    D_cello['C106']='C106[HUMANRECTALADENOCARCINOMA]'
    D_cello['CC20']='CC20[HUMANCOLONADENOCARCINOMA]'
    D_cello['NCIADRES']='NCIADRRES'
    D_cello['Hela']='HELA'
    D_cello['HL60(TB)']='HL60'
    
    keys=list(D_cello.keys())
    # Still missing
    for ccl_red in ccl_red:
        if not ccl_red in keys:
            print(f'{ccl_red}')
    #return([D_ccl[a] for a in ccls])
    
    return([D_cello[D_ccl[a]] for a in ccls])


def z_transform(M0):
    M0=M0-np.array(M0.mean(1))[:,None]
    M0=M0/np.array(M0.std(1))[:,None]
    return(M0)
def half_min_imputation(M0):
    mins=M0.min(axis=1)/2
    for i,j in zip(np.where(np.isnan(M0))[0], np.where(np.isnan(M0))[1]):
        M0.iloc[i,j]=mins[i]
    return (M0)


def load_atlantic_only_phos():
    # Phos+Full
    M0=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/lwn_imp_M0.pkl'))
    M0.columns=[a.split('_')[0]+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=sorted(set([a.split('_')[0] for a in D_prot.keys()]))
    
    return(D_prot, features, celllines)


def load_atlantic_only_full():
    # Phos+Full
    M0=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/als_imp_M1.pkl'))
    M0.columns=[a.split('_')[0]+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=sorted(set([a.split('_')[0] for a in D_prot.keys()]))
    
    return(D_prot, features, celllines)


def load_atlantic():
    # Phos+Full
    M0=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/lwn_imp_M0.pkl'))
    M1=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Atl/als_imp_M1.pkl'))
    M0=pd.concat([M0,M1]).copy()
    M0.columns=[a.split('_')[0]+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=sorted(set([a.split('_')[0] for a in D_prot.keys()]))
    
    return(D_prot, features, celllines)



def load_Gygi():
    M0=z_transform(pd.read_pickle('/home/icb/manuel.gander/Atl/n2/Imp/Gigy2/lwn_imp_M0.pkl'))
    M0.columns=[a+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=[a.split('_')[0] for a in D_prot.keys()]
    
    return(D_prot, features, celllines)

def load_kinase_scores():
    M0=pd.read_csv(f'{Path}/kinase_scores.tsv', sep='\t')
    M0.index=M0['PSP Kinases']
    M0=M0[M0.columns[128:-2]]
    M0=M0[[a for a in M0.columns if '_' in a]]
    M0.columns=map_celllines_to_cellosaurus([a.split('_')[0] for a in M0.columns])
    M0.columns=[a+f'_{i}' for i,a in enumerate(M0.columns)]
    fr=0.7
    M0=M0[np.isnan(M0).sum(1)<(1-fr)*len(M0.columns)].copy()
    M0=half_min_imputation(M0)
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=[a.split('_')[0] for a in D_prot.keys()]
    
    return(D_prot, features, celllines)

def load_RNA():
    M0=pd.read_pickle(f'{Path}/CCLE_RNA.pkl')
    M0.columns=[a+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=[a.split('_')[0] for a in D_prot.keys()]
    return(D_prot, features, celllines)

def load_identity():
    M0=pd.read_pickle(f'{Path}/Ccl_identity.pkl')
    M0.columns=[a+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=[a.split('_')[0] for a in D_prot.keys()]
    return(D_prot, features, celllines)

def load_drug_scores():
    M0=pd.read_csv(f'{Path}/drug_scores.tsv', sep='\t')
    M0.index=M0['Drug']
    M0=M0[M0.columns[128:-2]]
    M0=M0[[a for a in M0.columns if '_' in a]]
    M0.columns=map_celllines_to_cellosaurus([a.split('_')[0] for a in M0.columns])
    M0.columns=[a+f'_{i}' for i,a in enumerate(M0.columns)]
    D_prot={}
    for a in M0.columns:
        D_prot[a]=np.array(M0[a])
    features=list(M0.index)
    celllines=[a.split('_')[0] for a in D_prot.keys()]
    return(D_prot, features, celllines)

def load_dataset(dataset):
    if dataset=='atlantic':
        r = load_atlantic()
    elif dataset=='atl_only_phos':
        r = load_atlantic_only_phos()
    elif dataset=='atl_only_full':
        r = load_atlantic_only_full()
    elif dataset=='Gygi':
        r = load_Gygi()
    elif dataset=='kinase_scores':
        r = load_kinase_scores()
    elif dataset=='RNA':
        r = load_RNA()
    elif dataset=='drug_scores':
        r = load_drug_scores()
    elif dataset=='atl_full_and_kinases':
        r = load_atlantic_full_and_kinases()
    elif dataset=='identity':
        r = load_identity()
    else:
        print(f'Unknown dataset: {dataset}')
    a,b,c=r
    c=sorted(set(c))
    return a,b,c





def subsample_features(D_prot, features, celllines, fraction=0.5, seed=0):
    l=len(features)
    np.random.seed(seed)
    inds_kept = np.random.choice(list(range(l)), int(l*fraction), replace=False)
    for a in D_prot.keys():
        D_prot[a]=D_prot[a][inds_kept]
    features = np.array(features)[inds_kept]
    
    return(D_prot, features, celllines)

def subsample_celllines(D_prot, features, celllines, fraction=0.5, seed=0):
    ccls_i = sorted(set(D_prot.keys()))
    l=len(ccls_i)
    np.random.seed(seed)
    ccls_i_kept = np.random.choice(ccls_i, int(l*fraction), replace=False)
    D_prot2={}
    for a in ccls_i_kept:
        D_prot2[a]=D_prot[a]
    celllines=sorted(set(a.split('_')[0] for a in ccls_i_kept))
    
    return(D_prot2, features, celllines)


def prep_viability_AUCs2(D_prot):
    feature_data = list(D_prot.keys())
    dfv=pd.read_pickle(f"{Path}/Synched_aucs.pkl")
    dfv=dfv[np.isfinite(dfv['AUC'])].copy()
    frames=[]
    for cello in tqdm(sorted(set(dfv['Cello']))):
        cellos_found = [a for a in feature_data if a.split('_')[0]==cello]
        if len(cellos_found)>0:
            dfc=dfv[dfv['Cello']==cello].copy().copy()
            for c in cellos_found:
                dfc_copy=dfc.copy()
                dfc_copy['Cello_i']=c
                frames.append(dfc_copy)
    dfv=pd.concat(frames, ignore_index=True)
    import pickle
    # Load from pickle file
    with open(f'{Path}/Sensitivity_dict.pkl', 'rb') as file:
        Df = pickle.load(file)
    ls = list(dfv['Source'])
    ld = list(dfv['PubChem_CID'])
    lc = list(dfv['Cello'])
    
    senss = []
    for i in tqdm(range(len(dfv))):
        senss.append(Df[ls[i]].get(ld[i], {}).get(lc[i], 'nan'))
    dfv['sensitivity']=senss
    dfv = dfv[dfv['sensitivity']!='nan'].copy()
    return(dfv)


def get_holdouts_balanced(dfvs, reps, fr=0.2, source0='all'):
    if source0!='all':
        dfvs = dfvs[dfvs['Source']==source0].copy()
    all_ccls = sorted(set(dfvs['Cello']))
    sccls =  sorted(set(dfvs[dfvs['sensitivity']=='yes']['Cello']))
    nsccls = [a for a in all_ccls if not a in sccls]

    D_holdout = {}
    if len(sccls)>5:
        for i in range(reps):
            np.random.seed(i)
            holdouts_sens = np.random.choice(sccls, int(len(sccls)*fr))
            holdouts_non_sens = np.random.choice(nsccls, int(len(nsccls)*fr))
            D_holdout[i] = [*holdouts_non_sens, *holdouts_sens]
    return(D_holdout)

def do_gridpoint_for_dfvs(dfvs, D_prot, kn, alpha, l1_ratio, reps, source0='all'):
    ccls=sorted(set(dfvs['Cello']))
    res={}
    if len(ccls)>30:
        D_holdout = get_holdouts_balanced(dfvs, reps, source0=source0)
        if len(D_holdout)>0:
            res['pearsons']=[]
            res['spearmans']=[]
            res['l1_ratios']=[]
            res['RMSE']=[]
            res['RMSE_mean_model']=[]
            for hs in range(reps):
                holdouts=D_holdout[hs]
                dfv0, dfh=split_of_validation_cell_lines(dfvs, ccls=holdouts)
                train_input=dfv_to_train_arrays(dfv0, D_prot)
                mean = train_input[1].mean()
        
                # Preselection of kn features and inference
                predictions, ys = calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio)
                p=predictions[0]
                y=ys[0]
        
                res['pearsons'].append(scipy.stats.pearsonr(p, y).statistic)
                res['spearmans'].append(scipy.stats.spearmanr(p, y).statistic)
                res['RMSE'].append(np.sqrt(np.mean((p-y)**2)))
                res['l1_ratios'].append(np.mean(abs(p-y))/np.mean(abs(mean-y)))
                res['RMSE_mean_model'].append(np.sqrt(np.mean((mean-y)**2)))
    return(res)


def calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio):
    train_input=dfv_to_train_arrays(dfv0, D_prot)    
    rss1 = paerson_nb(train_input[0], train_input[1])

    hold_input=dfv_to_train_arrays(dfh, D_prot)
    
    
    if kn<len(rss1):
        r_kept=np.argsort(rss1)[-kn:]
        X_kept=train_input[0][:,r_kept]
        Xh_kept=hold_input[0][:,r_kept]
    else:
        X_kept=train_input[0]
        Xh_kept=hold_input[0]
    y_train=train_input[1]
    y_test = hold_input[1]
                                  
    regr = ElasticNet(random_state=0, alpha=alpha, l1_ratio=l1_ratio)
    regr.fit(X_kept, y_train)
    
    prediction = regr.predict(Xh_kept)
    
    return([prediction], [y_test])



def do_grid_point(drug, kn, l1_ratio, alpha, dfv, D_prot, dataset, reps=100):
    #try:
    dfg=dfv[dfv['PubChem_CID']==drug].copy()
    ccls=sorted(set(dfg['Cello']))
    if len(ccls)>30:
        # create some holdouts
        D_holdout={}
        for i in range(reps):
            np.random.seed(i)
            D_holdout[i]=np.random.choice(ccls, int(len(set(dfg['Cello']))*0.15))

        res={}
        res['pearsons']=[]
        res['spearmans']=[]
        res['l1_ratios']=[]
        for hs in range(reps):
            holdouts=D_holdout[hs]
            dfv0, dfh=split_of_validation_cell_lines(dfg, ccls=holdouts)
            train_input=dfv_to_train_arrays(dfv0, D_prot)
            mean = train_input[1].mean()

            # Preselection of kn features and inference
            predictions, ys = calc_predictions_combinations(dfv0, dfh, D_prot, kn, alpha, l1_ratio)
            p=predictions[0]
            y=ys[0]

            res['pearsons'].append(scipy.stats.pearsonr(p, y).statistic)
            res['spearmans'].append(scipy.stats.spearmanr(p, y).statistic)
            res['l1_ratios'].append(np.mean(abs(p-y))/np.mean(abs(mean-y)))
    else:   
        res={}
        res['pearsons']=[np.NaN]
        res['spearmans']=[np.NaN]
        res['l1_ratios']=[np.NaN]
    #except:
     #   res={}
     #   res['pearsons']=[np.NaN]
     #   res['spearmans']=[np.NaN]
     #   res['l1_ratios']=[np.NaN]
        
    dfr=pd.DataFrame({'Drug':drug, 'kn':kn, 'l1_ratio':l1_ratio, 'alpha':alpha, 'dataset':dataset, 'n_ccls':len(ccls),
                     'pearsons':res['pearsons'], 'spearmans':res['spearmans'], 'result_l1_ratios':res['l1_ratios'],
                     'holdout_set':range(len(res['spearmans']))})
    return dfr
