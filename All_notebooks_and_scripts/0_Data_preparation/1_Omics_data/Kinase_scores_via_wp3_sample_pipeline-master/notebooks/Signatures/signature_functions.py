
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

#entity_subtypes = 'tissue_topology'
entity_subtypes = 'Histologic Subtype'



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


def get_final_signatures_from_t_test_results_by_criteria(signatures_df: pd.DataFrame,
                                                         p_value_threshold: float = 0.01,
                                                         fdr_threshold: float = .01,
                                                         delta_threshold: float = .75) -> pd.DataFrame:
    """
    Returns the final list of the signatures after filtering based on the selected criteria
    """
    signatures_df['delta'] = abs(signatures_df['means_group1'] - signatures_df['means_group2'])
    if p_value_threshold:
        signatures_df = signatures_df[(signatures_df['p_values'] < p_value_threshold)]

    if fdr_threshold:
        signatures_df = signatures_df[(signatures_df['fdr'] < fdr_threshold)]

    if delta_threshold:
        signatures_df = signatures_df[signatures_df['delta'] > delta_threshold]

    return signatures_df


def intersection(lst1: list, lst2: list):
    """
    The intersection between two lists
    """
    return list(set(lst1) & set(lst2))


def set_diff(lst1: list, lst2: list):
    """
     Items in list1 which do not exist in list2
    """
    return list(set(lst1) - set(lst2))


def read_intensity_file(path_to_file: str,index_col:str) -> pd.DataFrame:
    """ Report directory is the output of the TOPAS WP3 pipeline"""
    intensity_file = pd.read_csv(path_to_file)
    intensity_file = intensity_file.set_index(index_col)
    return intensity_file


def subset_intensity_file_fp(report_directory: str,
                             meta_data: pd.DataFrame,
                             selectedProteins: list,
                             selectedPatientsEntity: str) -> pd.DataFrame:
    """
    Subsets the intensity file based on selected entity and selected proteins
    """
    path_to_file = f'{report_directory}/preprocessed_fp.csv'
    intensity_file = pd.read_csv(path_to_file)

    intensity_file = intensity_file[intensity_file['Gene names'].isin(selectedProteins)]
    selected_patients = list(meta_data.loc[meta_data[entity_subtypes] == selectedPatientsEntity, 'Sample name'])
    selected_df = intensity_file[intersection(intensity_file.columns, selected_patients)]
    selected_df.index = selectedProteins
    return selected_df


def prepareDataframeforTest(intensity_file:pd.DataFrame,
                            meta: pd.DataFrame,
                            minimum_patients_per_entity: int = 20,
                            protein_expressed_in_at_least_percent: int = 70):
    """
    INPUT: - report_directory: is the result folder of the wp3 clinical proteomics pipeline
           - fp_pp: fp for full protein and pp for phspho proteins
           - metapath: is the path to the meta file
           - minimum_patients_per_entity: the minimum patients in the entity to consider it as entity
           - protein_expressed_in_at_least_percent: the protein expressed within at least this percentage of all patients
           - index_col: the column in the intensity to select as the index of the intensity file

    OUTPUTs: - a dataframe before T_test or ANOVA
             - a list of the protein names

    """

    entity_count = meta.groupby(entity_subtypes)['Sample name'].count()
    entities = entity_count[entity_count >= minimum_patients_per_entity]  # only the entities with more than n patienst will be used
    entities = list(entities.index)
    meta_data = meta[meta[entity_subtypes].isin(entities)]
    meta_data = meta_data[['Sample name', entity_subtypes]]

    # Filtering and merging the intensity table to that of the meta data
    patients_for_ANOVA = intersection(intensity_file.columns, meta_data['Sample name'])
    final_inensity = intensity_file[patients_for_ANOVA].transpose()
    final_inensity.columns = intensity_file.index
    # final_inensity = final_inensity.dropna(how='all',axis=1)
    final_inensity['Sample name'] = final_inensity.index
    final_df = pd.merge(final_inensity, meta_data, on='Sample name', validate='one_to_one')
    protein_peptide = list(final_inensity.columns)
    protein_peptide.remove('Sample name')
    subset = final_df[protein_peptide]
    protein_peptide = subset.columns[subset.count(axis=0) > ((protein_expressed_in_at_least_percent * 0.01) * len(subset.index))].to_list()
    return final_df, protein_peptide


def anova_test(inputDF: pd.DataFrame,
               protein_peptide: list,
               metaDataColumn: str) -> pd.DataFrame:
    """
    Performs ANOVA for each protein/peptide row wise
    INPUT: - a dataframe where the protein/peptides are the rows and
                               the patients are the columns
           - protein_peptide is the list of the columns to do the ANOVA
           - metaDatacolumns is the column in the inputDf which contains the subtypes for grouping

    """
    # new_Df = new_Df.transpose()
    F_tests, p_values, mean_grps = [], [], []
    for gene in protein_peptide:
        f, p = 1, 1  # a precomputed value for the p and F
        column_names = [gene, metaDataColumn]
        grp = inputDF[column_names]
        grp = grp.dropna()
        #
        result = grp.groupby(metaDataColumn)[gene].apply(list)
        average = grp.groupby(metaDataColumn)[gene].mean()
        average = list(average)
        # more than two groups for each gene
        if len(grp['localisation'].unique()) > 2:
            f, p = f_oneway(*result)
        F_tests.append(f)
        p_values.append(p)
        mean_grps.append(average)
    p_df = pd.DataFrame(list(zip(F_tests, p_values, mean_grps)))
    p_df.columns = ['F_tests', 'p_values', 'means_pergroup']
    p_df['Gene Names'] = protein_peptide
    p_df = p_df[p_df['p_values'].notna()]
    fdr_multi_correction = fdrcorrection(p_df.p_values, alpha=0.05, method='indep', is_sorted=False)
    p_df['fdr'] = list(fdr_multi_correction[1])
    return p_df


def t_test(x: tuple, y: tuple):
    """ Performs t_test from two numerical tuples as the two groups"""
    t_test = ttest_ind(x, y)
    F = list(t_test)[0]
    p = list(t_test)[1]
    return F, p


def one_vs_all_t_test(inputDF: pd.DataFrame,
                      protein_peptide: list,
                      favoriteentity: str,
                      metaDataColumn: str) -> pd.DataFrame:
    """
    Performs a one vs all T_test per each protein/peptide for the patients with the favorite entity vs all other entities
    to find out the signature protein for that specific entity

    :inputDF: a dataframe where the patients are the rows and
                               the proteins/peptides are the columns
    :protein_peptide: the list of the columns to do the t_test based on
    :favoriteentity: the main group for the t_test i.e: chordoma
    :metaDatacolumns: the column in the inputDf which contains the subtypes for grouping

    """
    F_tests, p_values, mean_grp1, mean_grp2, num_grp1, num_other_grps = [], [], [], [], [], []
    for gene in protein_peptide:
        f, p = 1, 1  # a precomputed value for the p and F
        column_names = [gene, metaDataColumn]
        grp = inputDF[column_names]
        grp = grp.dropna()
        df = grp
        g1 = df[(df[metaDataColumn] == favoriteentity)]  # first group for the test
        g2 = df[(df[metaDataColumn] != favoriteentity)]  # second group for the test
        x = tuple(g1[gene])
        y = tuple(g2[gene])
        f, p = t_test(x, y)
        F_tests.append(f)
        p_values.append(p)
        mean_x = g1[gene].mean()
        mean_y = g2[gene].mean()
        num_grp1.append(len(g1[gene]))
        num_other_grps.append(len(g2[gene]))
        mean_grp1.append(mean_x)
        mean_grp2.append(mean_y)
    p_df = pd.DataFrame(list(zip(F_tests, p_values, mean_grp1, mean_grp2, num_grp1, num_other_grps)))
    p_df.columns = ['t_statistics', 'p_values', 'means_group1', 'means_group2', 'num_samples_groups_interest', 'num_sample_other_groups']
    p_df['Gene Names'] = protein_peptide
    p_df = p_df[p_df['p_values'].notna()]
    fdr_multi_correction = fdrcorrection(p_df.p_values, alpha=0.05, method='indep', is_sorted=False)
    p_df['fdr'] = list(fdr_multi_correction[1])
    p_df['up_down'] = 'up'
    p_df['up_down'][(p_df['means_group1'] < p_df['means_group2'])] = 'down'
    return p_df


def get_weights_per_protein_using_average_entity_intensities(intensity_file: pd.DataFrame,
                                                             meta_data: pd.DataFrame,
                                                             list_signature_proteins: pd.Series,
                                                             entity: str,
                                                             list_patients: pd.Series):
    """
    Returns the vector of weights  based on the average of intensities for the patients of interested  cohort
    compared to other patients it can  be either +1 or -1.


    """
    intensity_file = intensity_file.loc[list_signature_proteins, :]
    # calculation of the weights
    positive_pateints = get_list_positive_labeled_patients(meta_data, list_patients, entity)
    weights = [0 for i in range(len(intensity_file.index))]
    for j in range(len(weights)):
        protein = intensity_file.index[j]
        grp = intensity_file[intensity_file.index == protein]
        grp_mean1 = grp[positive_pateints].mean(axis=1, skipna=True)  # mean of positive patints
        grp_mean2 = grp.drop(positive_pateints, axis=1).mean(axis=1, skipna=True)  # mean of the negative patients
        if grp_mean1[0] > grp_mean2[0]:
            weights[j] = 1  # 1 for up regulation
        else:
            weights[j] = -1  # -1 for down regulation
    return weights


def get_list_positive_labeled_patients(meta_data: pd.DataFrame,
                                       patients_list: list,
                                       entity: str):
    positive_labeled_patients = list(meta_data.loc[meta_data[entity_subtypes] == entity, 'Sample name'])
    positive_labeled_patients = intersection(positive_labeled_patients, patients_list)
    return positive_labeled_patients




def cross_validation_for_one_vs_all_t_test(inputDF: pd.DataFrame,
                                           intensity_file: pd.DataFrame,
                                           meta_data: pd.DataFrame,
                                           entity: str,
                                           protein_peptides_list,
                                           k_folds: int = 5,
                                           num_repeats: int = 1,
                                           p_value_cutoff: float = .01,
                                           fdr_cutoff: float = .01,
                                           average_difference_twogroups: float = .75) -> pd.DataFrame:

    """
    inputDF: the inputDF which is also used for the t_test  and is the output of prepareDataframeforTest
    raw_df: preprocessed_fp.csv from the wp3 pipeline
    """
    x = inputDF.to_numpy()
    y = inputDF[entity_subtypes].to_numpy()

    skf = RepeatedStratifiedKFold(n_splits=k_folds, n_repeats=num_repeats)
    all_data = []
    for train_index, test_index in skf.split(x, y):
        # training: calculation of the signatures based on one_vs_all_t_test
        train_set = inputDF.iloc[train_index, :]
        fp_signatures = one_vs_all_t_test(train_set, protein_peptides_list, entity, entity_subtypes)
        primary_signatures = get_final_signatures_from_t_test_results_by_criteria(fp_signatures, p_value_cutoff, fdr_cutoff,
                                                                                  average_difference_twogroups)
        primary_signatures_proteins = list(primary_signatures['Gene Names'])

        # test: calculation of the disease_score
        test_set = inputDF.iloc[test_index, :]
        test_patients = intersection(intensity_file.columns, test_set['Sample name'])
        sub_intensity = intensity_file[intensity_file.index.isin(primary_signatures_proteins)]
        sub_intensity = sub_intensity[test_patients]

        weights_array = get_weights_per_protein_using_average_entity_intensities(sub_intensity, meta_data, primary_signatures_proteins,
                                                                                 entity, test_patients)
        positive_labeled_patients = get_list_positive_labeled_patients(meta_data, test_patients, entity)
        results = calculate_entity_score_from_the_signature_proteins(sub_intensity, positive_labeled_patients, weights_array,
                                                                     transpose_weights=True)
        all_data.append(results)

    merged_df = pd.concat(all_data)
    return merged_df


def univariate_ROC_analysis_by_CV_permutation(pre: pd.DataFrame,
                                              entity: str,
                                              kFold: int = 5,
                                              repeats: int = 10,
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
    skf = RepeatedStratifiedKFold(n_splits=kFold, n_repeats=repeats)
    aucs = []
    for train_index, test_index in skf.split(x, y):
        auc = ROC_curve_analysis(y[train_index], x[train_index], curve_title='', plot=False)
        aucs.append(auc)
    aucs = pd.DataFrame(aucs)
    aucs = aucs.dropna()
    aucs.columns = ['auc']
    fulfilled = len(aucs[(aucs['auc'] > threshold) | (aucs['auc'] < (1 - threshold))])
    return (fulfilled / (kFold * repeats)) * 100



def univariate_analysis_with_CV_for_all_proteins(inputDF:pd.DataFrame,
                                                list_proteins:pd.DataFrame,
                                                entity,
                                                K_folds:int = 5,
                                                repeats:int = 1,
                                                threshold:float = 0.7
                                                ) -> pd.DataFrame:
    """
    Performs univariate CV analysis 
    :inputDF: the resutl of the preapre function
    :list_proteins : the list of the proteins to perform the test ob
    :entity: the positive labes entity
    """


    roc_res = [0 for i in range(len(list_proteins))]
    for i in range(len(list_proteins)):
        df = inputDF.loc[:, [list_proteins[i], entity_subtypes]].dropna()
        if (len(df.index) > K_folds):
            try:
                roc_res[i] = univariate_ROC_analysis_by_CV_permutation(df, entity, kFold=K_folds,
                                                                            repeats=repeats, threshold=threshold,
                                                                            scores=str(list_proteins[i]),
                                                                            labels=entity_subtypes)
            except:
                print(f'Error for the {list_proteins[i]}')
    
    return pd.DataFrame(list(zip(list_proteins, roc_res)),columns=['names',f'CV_robustness_above_{threshold}'])





def calculate_entity_score_from_the_signature_proteins(intensity_file: pd.DataFrame,
                                                       positive_labeled_patients: list,
                                                       weights_list: list,
                                                       transpose_weights: bool = False) -> pd.DataFrame:
    """
    The prediction is based on the model:     scores = [intensity] * W
    W is the vector of weights
    cohort with respect to other patients it can  be either +1 or -1.
    The intensities NAs will be imputed based on the average of patients for that protein

    :intensity_file: A dataframe with rows as protein and the columns as the patients
    :weights_array: list of the labels
    :real_labels: list of the 1 for positive and 0 for the negative labeled patients
    OUTPUT:
        dataframe of
        - labels a vector of 0 and 1, real values if the patient is with the cohort it will be 1 other wise it will be 0
        - scores the predicted result of [intensity] * weights
    """
    labels = [1 if i in positive_labeled_patients else 0 for i in intensity_file.columns]
    intensity_file = intensity_file.fillna(
        intensity_file.mean())  # imputation for each gene as the average of the intensities across all patients
    weights = np.array(weights_list)
    if transpose_weights:
        weights = weights.transpose()
    scores = np.dot(intensity_file.transpose().to_numpy(), weights)
    result_df = pd.DataFrame(list(zip(labels, scores)))
    result_df.columns = ['labels', 'scores']
    result_df.index = intensity_file.columns
    return result_df


def hyper_parameters_tuning_by_CV_one_vs_all_SVM(raw_df: pd.DataFrame,
                                                protein_peptides_list: list,
                                                entity: str) -> dict:
    """
    Returns the best tuned Hyper parameters as a dictionary of C and kernel: 'linear'
    for the one vs all classification task using support vector machines by Grid Search
    it use the linear kernel function to get the coefficints
    : raw_df: Data frame with the patients as rows and the proteins as the columns an exta columns should be named Sarcomasubtype with the class labels
    : protein_peptides_list: is the list of the proteins for the tuning
    : Entity: is the name of the Entity such as Chordoma, Sarcoma etc
    """
    input_svm = raw_df
    y = input_svm[entity_subtypes].to_numpy()
    y[y != entity] = f'not_{entity}'  # negative class; so we only have Entity and not_Entity
    input_svm = input_svm[protein_peptides_list]
    imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    x = imputer.fit_transform(input_svm)  # column-wise average imputation
    svc = SVC()
    # we only use the linear kernel to use it as a feature selection method
    parameters = [
        {'C': [1, 2, 5, 10, 100, 1000], 'kernel': ['linear']}
        # {'C':[1, 10, 100, 1000], 'kernel':['rbf'], 'gamma':[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]},
        # {'C':[1, 10, 100, 1000], 'kernel':['poly'], 'degree': [2,3,4] ,'gamma':[0.01,0.02,0.03,0.04,0.05]}
    ]
    grid_search = GridSearchCV(estimator=svc,
                               param_grid=parameters,
                               scoring='accuracy',
                               cv=5,
                               verbose=0)
    model = grid_search.fit(x, y)
    return grid_search.best_params_


def predict_by_SVM_one_vs_all(inputSVM: pd.DataFrame,
                              entity: str,
                              labels_column: str,
                              model_params: dict,
                              protein_peptides_list: list,
                              impute_strategy='mean') -> dict:
    """
    Returns a dictionary of the coefficients, predicted classed and confusion matrix after training classification task using SVM
    :inputSVM: data frame with the rows as patients and the columns as the proteins + class labels as a separate column
    :Entity: the positive label class like Chordoma
    :labels_column: the columns with two class labels (Entity is mandatory and non-entity(any name could be possible))
    :model_params: linear kernel + C hyper param as a dictionary
    :protein_peptides_list:list of the proteins to do the classification
    """
    y = inputSVM[labels_column].to_numpy()
    y[y != entity] = f'not_{entity}'  # negative class; so we only have Entity and not_Entity
    inputSVM = inputSVM[protein_peptides_list]
    imputer = SimpleImputer(missing_values=np.nan, strategy=impute_strategy)
    x = imputer.fit_transform(inputSVM)  # column-wise average imputation
    trained_svm = SVC(**model_params).fit(x, y)
    prediction_svm = trained_svm.predict(x)
    confusionmatrix = confusion_matrix(y, prediction_svm)
    classificationreport = classification_report(y, prediction_svm)
    f1score = f1_score(y, prediction_svm, average='macro')
    coefs = trained_svm.coef_
    return {'coefficients': coefs,
            'predicted_labels': prediction_svm,
            'confusion_matrix': confusionmatrix,
            'f1_score': f1score,
            'classification_report': classificationreport}



def SVM_cross_validaition_with_t_test_Feature_selection(input_df: pd.DataFrame,
                                                        intensity_file:pd.DataFrame,
                                                        meta_df: pd.DataFrame,
                                                        list_proteins:list,
                                                        Entity: str,
                                                        p_value_cutoff = None,
                                                        fdr_cutoff = .01,  # for t_test as the feature selecetor
                                                        average_difference_twogroups = .75,  # the difference between the two groups
                                                        n_splits = 5,
                                                        n_repeats = 1,
                                                        weight_multiplier = -1
                                                        ) -> pd.DataFrame:
    """
    Splits the data into k-folds
    for each fold the train set will be the intensities/patwayscores
    the response will be a binary vector with Entity as the + class
    and all other entities as the second class
    for each fold performs t_test to select the signatures uses the selected features 
    for parameter hyper tuning with SVM and uses the selected paramters and feeatures for the 
    testset to predict entity scores
    :input_df: Dataframe with rows as patients and columns as features two extra column "Sample name" metadata 
    :intensity_file: an indexed dataframe with patients as index and the features are the columns
    :list_proteins: list of the proteins/peptides/kinases/baskets
    """
    # selecting the features based on the 
    Y = input_df[entity_subtypes].to_numpy()
    X = input_df.to_numpy()
    Entity_patients = get_list_positive_labeled_patients(meta_df, input_df['Sample name'], Entity)
    skf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats)
    all_data = []
    for train_index, test_index in skf.split(X, Y):
        train_set = input_df.iloc[train_index, :]

        # FEATURE SELECTION: by t_test with train set
        fp_signatures = one_vs_all_t_test(train_set,
                                            list_proteins,
                                            Entity,
                                            entity_subtypes)
        
        primary_signatures = get_final_signatures_from_t_test_results_by_criteria(fp_signatures,
                                                                                    p_value_cutoff,
                                                                                    fdr_cutoff,
                                                                                    average_difference_twogroups)
        primary_signatures_proteins = list(primary_signatures['Gene Names'])

        # training by support vector machines
        best_params = hyper_parameters_tuning_by_CV_one_vs_all_SVM(train_set,
                                                                    primary_signatures_proteins,
                                                                    Entity)
        
        trained_model = predict_by_SVM_one_vs_all(train_set,
                                                    Entity,
                                                    entity_subtypes,
                                                    best_params,
                                                    primary_signatures_proteins)

        # testset with the selected features
        test_set = input_df.iloc[test_index, :]
        sub_intensity = intensity_file[intensity_file.index.isin(primary_signatures_proteins)]
        sub_intensity = sub_intensity[test_set['Sample name']]
        weights = trained_model['coefficients'] * weight_multiplier
        scores = calculate_entity_score_from_the_signature_proteins(sub_intensity,
                                                                    Entity_patients,
                                                                    weights,
                                                                    transpose_weights=True)
        all_data.append(scores)
    all_data = pd.concat(all_data)
    ROC_curve_analysis(all_data['labels'], all_data['scores'], '')
    return all_data


def config_reader(config_path: str) ->dict:
    with open(config_path, 'r') as f:
        config = json.load(f)
    return config





