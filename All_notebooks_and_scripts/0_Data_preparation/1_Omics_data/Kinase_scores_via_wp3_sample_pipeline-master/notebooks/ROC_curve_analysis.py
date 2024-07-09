
import pandas as pd
import json
import numpy as np

from scipy.stats import f_oneway, ttest_ind
from statsmodels.stats.multitest import fdrcorrection

from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV
from sklearn import metrics
from sklearn.svm import SVC
from sklearn.impute import SimpleImputer
from sklearn.metrics import confusion_matrix, classification_report, f1_score

import matplotlib.pyplot as plt



def ROC_curve_analysis(labels: np.ndarray,
                       scores: np.ndarray,
                       curve_title: str,
                       plot: bool = True):
    """
    Plots the ROC curve using matplotlib
    labels: 1 and 0 indicating positive and negative classes
    scores: continuous variable the more value is favored for positive class
    """
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores)
    roc_auc = metrics.auc(fpr, tpr)
    if plot:
        display = metrics.RocCurveDisplay(fpr=fpr,
                                          tpr=tpr,
                                          roc_auc=roc_auc,
                                          estimator_name=curve_title)
        display.plot()
        plt.show()
    return roc_auc


def univariate_ROC_analysis_by_CV_permutation(pre: pd.DataFrame,
                                              entity: str,
                                              kFold: int = 5,
                                              repeats: int = 20,
                                              threshold: float = 0.5,
                                              scores: str = 'scores',
                                              labels: str = 'labels') -> float:
    """
    Return the percentage of the stability for the AUCs of ROC above the threshold value
    :pre: A dataframe with the two columns scores as numerical and labels as strings
    :entity: positive class labels
    """
    x = pre[scores].copy().to_numpy()
    y = pre[labels].copy()
    y[y == entity] = 1
    y[y != 1] = 0
    y = y.astype(int).to_numpy()
    ROC_curve_analysis(y, x, curve_title='ROC analysis_Sensitive VS. Resistant', plot=True)
    skf = RepeatedStratifiedKFold(n_splits=kFold, n_repeats=repeats)
    aucs = []
    for train_index, test_index in skf.split(x, y):
        auc = ROC_curve_analysis(y[train_index], x[train_index], curve_title='', plot=False)
        aucs.append(auc)
    aucs = pd.DataFrame(aucs)
    aucs = aucs.dropna()
    aucs.columns = ['auc']
    fulfilled = len(aucs[(aucs['auc'] > threshold) | (aucs['auc'] < (1 - threshold))])
    return aucs

df = pd.read_excel('/home/amir/Desktop/Annika_files/TUPAC score_ROC curve.xlsx')
Aucs_df = univariate_ROC_analysis_by_CV_permutation(df,'S',scores='EGFR',labels='Response')
Aucs_df.to_excel('/home/amir/Desktop/Annika_files/ROC_EGFR_R_S.xlsx',index=False)