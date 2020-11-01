#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:22:25 2020

@author: shanjuyeh
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import preprocessing
from sklearn.model_selection import RandomizedSearchCV
import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, roc_auc_score, precision_score, recall_score, f1_score, \
    confusion_matrix, matthews_corrcoef, precision_recall_fscore_support, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from keras.utils import to_categorical
from keras.layers import Input, Dense, Dropout
from keras.models import Model
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from keras import regularizers
from keras import Sequential
from keras.callbacks import EarlyStopping
from keras import backend as K
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from keras.wrappers.scikit_learn import KerasClassifier
import matplotlib.pyplot as plt
import random
import imblearn
from imblearn.over_sampling import SMOTE
import pickle

meta = pd.read_csv('./input/age_meta_GPL570.csv')
data = pd.read_csv('./input/GPL570_final.csv')
data = data.fillna(0)
data = data.rename(columns={'Unnamed: 0': 'gsm'})
data = data.drop(['DDX3Y', 'EIF1AY', 'KDM5D', 'RPS4Y1', 'USP9Y'], axis=1) # remove gender features

data_used = data[ (data['age'] == 0) |  (data['age'] == 1) | (data['age'] == 2)]
data_used = data_used.set_index('id')
data_used = data_used.drop(['tissue', 'gender'], axis=1)

data_used = data_used.drop('age', axis=1)
used_id = pd.DataFrame(data_used.index)
used_id = used_id.rename(columns={'id':'gsm'})

meta_used = pd.merge(used_id, meta, on='gsm')

matrix = data_used
label = meta_used['label']

x = matrix.to_numpy()
y = label.to_numpy()

# Indicies of each class observations
i_0 = np.where(y == 0)[0]
i_1 = np.where(y == 1)[0]
i_2 = np.where(y == 2)[0]

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.33)

# Oversampling by SMOTE
oversample = SMOTE()
x_trn, y_trn = oversample.fit_resample(x_train, y_train)

# Downsampleing
i_class0 = np.where(y_trn == 0)[0]
i_class1 = np.where(y_trn == 1)[0]
i_class2 = np.where(y_trn == 2)[0]
n_class0 = 6500
n_class1 = 7000
n_class2 = 6500
i_class0_downsampled = np.random.choice(i_class0, size=n_class0, replace=False)
i_class1_downsampled = np.random.choice(i_class1, size=n_class1, replace=False)
i_class2_downsampled = np.random.choice(i_class2, size=n_class2, replace=False)

yy_trn = np.hstack((y_trn[i_class1_downsampled], y_trn[i_class2_downsampled], y_trn[i_class0_downsampled]))
xx_trn = np.vstack((x_trn[i_class1_downsampled], x_trn[i_class2_downsampled], x_trn[i_class0_downsampled]))


# XGBOOST
params = {
        'min_child_weight': [1, 5, 10],
        'gamma': [0.5, 1, 1.5, 2, 5],
        'subsample': [0.6, 0.8, 1.0],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'max_depth': [5, 10, 20]
        }

XGB = XGBClassifier(learning_rate=0.03, n_estimators=1300, silent=True, nthread=1)
folds = 5
param_comb = 5
skf = StratifiedKFold(n_splits=folds, shuffle = True, random_state = 1001)
model = RandomizedSearchCV(XGB, param_distributions=params, n_iter=param_comb, n_jobs=4, cv=skf.split(xx_trn, yy_trn), verbose=3, random_state=1001 )
model.fit(xx_trn, yy_trn)

y_pred = model.predict(x_test)
target_names = ['group1', 'group2', 'group3']
report = classification_report(y_test, y_pred, target_names=target_names)

std = np.std(model.cv_results_['mean_test_score'])
print("Accuracy of model is: ", accuracy_score(y_test, y_pred))
print("Precision = {}".format(precision_score(y_test, y_pred, average='weighted')))
print("Recall = {}".format(recall_score(y_test, y_pred, average='weighted')))

## Prediction Score ============================================
un_GPL = data[data['age']== -1]
ID = un_GPL[['id']]
ID['age_old']='NA'
ID = ID.reset_index(drop=True)
ID = ID.rename(columns={'age_old':'age'})

un_GPL_used = un_GPL.drop(['age', 'tissue', 'gender'], axis=1)
un_GPL_used = un_GPL_used.set_index('id')
un_GPL_used = un_GPL_used.fillna(0)
un_x = un_GPL_used.to_numpy()

y_score_un = model.predict_proba(un_x)
y_pred_un = model.predict(un_x)
y_pred_un = pd.DataFrame(y_pred_un)
y_pred_un = y_pred_un.rename(columns={0:'pred'})

y_score_un = pd.DataFrame(y_score_un)
y_score_un = y_score_un.rename(columns={0:'group1', 1:'group2', 2:'group3'})
GPL_un_finish = pd.concat([ID, y_pred_un, y_score_un], axis=1)

y_score = model.predict_proba(x)
y_pred = model.predict(x)
y_pred = pd.DataFrame(y_pred)
y_pred = y_pred.rename(columns={0:'pred'})
y_score = pd.DataFrame(y_score)
y_score = y_score.rename(columns={0:'group1', 1:'group2', 2:'group3'})
com = meta_used[['gsm', 'label']]
GPL_finish = pd.concat([com, y_pred, y_score], axis=1)
GPL_finish = GPL_finish.rename(columns={'gsm':'id', 'label':'age'})

GPL_complete = pd.concat([GPL_finish, GPL_un_finish], axis=0)
GPL_complete.to_csv('./output/GPL570_result_1.csv')

g1 = GPL_complete[GPL_complete['pred']==0].shape[0]
g2 = GPL_complete[GPL_complete['pred']==1].shape[0]
g3 = GPL_complete[GPL_complete['pred']==2].shape[0]

a = g1/142257
b = g2/142257
c = g3/142257

print(a)
print(b)
print(c)







