import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
import pickle as pkl
import os


def read_data(dataset):
    fn = '../output/LABELS/' + dataset + '_final.csv'
    df = pd.read_csv(fn)
    return df

def split_data(df, dataset, n_task=2):
    gender = df['gender'].to_numpy(dtype=np.int32)
    tissue = df['tissue'].to_numpy(dtype=np.int32)
    age = df['age'].to_numpy(dtype=np.int32)

    feature = df.iloc[:, 4:].to_numpy(dtype=np.float32)
    print(type(gender), feature.shape, np.unique(tissue, return_counts=True))

    if n_task==2:
        # [index gender, index tissue]
        idx = np.concatenate([np.array(gender != -1).reshape(-1, 1), np.array(tissue != -1).reshape(-1, 1)], axis=1)
        print(idx)
    else:
        # [index gender, index tissue, index age]
        idx = np.concatenate([np.array(gender != -1).reshape(-1, 1), np.array(tissue != -1).reshape(-1, 1),
                              np.array(age != -1).reshape(-1, 1)], axis=1)
        print(idx)

    labeled_any_idx = np.any(idx, axis=1)
    labeled_all_idx = np.all(idx, axis=1)
    print(sum(labeled_any_idx), sum(labeled_all_idx))

    feature_labeled_any = feature[labeled_any_idx, :]
    gender_labeled_any = gender[labeled_any_idx]
    tissue_labeled_any = tissue[labeled_any_idx]
    age_labeled_any = age[labeled_any_idx]

    skf = StratifiedKFold(n_splits=5, random_state=0)
    fold = 0
    os.makedirs('../output/LABELS/' + dataset + '_task_' + str(n_task), exist_ok=True)
    for train_index, test_index in skf.split(feature_labeled_any, tissue_labeled_any):
        x_train, x_test = feature_labeled_any[train_index, :], feature_labeled_any[test_index, :]
        y_train_g, y_train_t, y_train_a, \
        y_test_g, y_test_t, y_test_a = gender_labeled_any[train_index], tissue_labeled_any[train_index], age_labeled_any[train_index],\
                                       gender_labeled_any[test_index], tissue_labeled_any[test_index], age_labeled_any[test_index]

        fn = '../output/LABELS/' + dataset + '_task_' + str(n_task) + '/fold_' + str(fold) + '.pkl'
        with open(fn, 'wb') as f:
            if n_task==2:
                pkl.dump([x_train, y_train_g, y_train_t, x_test, y_test_g, y_test_t], f)
            else:
                pkl.dump([x_train, y_train_g, y_train_t, y_train_a, x_test, y_test_g, y_test_t, y_test_a], f)

        fold += 1

    fn = '../output/LABELS/' + dataset + '_task_' + str(n_task) + '/fold_all.pkl'
    with open(fn, 'wb') as f:
        if n_task==2:
            pkl.dump([feature_labeled_any, gender_labeled_any, tissue_labeled_any, feature, gender, tissue], f)
        else:
            pkl.dump([feature_labeled_any, gender_labeled_any, tissue_labeled_any, age_labeled_any, feature, gender, tissue, age], f)


def cv_result(dataset, n_task):
    path = '../output/LABELS/' + dataset + '_task_' + str(n_task)
    dropout_list = [0.1, 0.2, 0.3, 0.4, 0.5]
    accs = []
    stds = []
    for i in range(len(dropout_list)):
        do = dropout_list[i]
        v = []
        for j in range(5):
            fn = path + '/dropout_' + str(do) + '_fold_' + str(j) + '.txt'
            with open(fn, 'r') as f:
                for last in f:
                    pass
                line = last.strip().split(' ')
                print(line)
            acc = np.array([float(line[n_task+1]), float(line[n_task+2])])
            v.append(acc)
        v = np.array(v)
        v_avg = np.mean(v, axis=0)
        std_avg = np.std(v, axis=0)
        accs.append(v_avg)
        stds.append(std_avg)

    accs = np.array(accs)
    stds = np.array(stds)
    idx = np.argmax(np.sum(accs, axis=1))
    print(accs)
    print(stds)
    print('best dropout rate: {:.1f}'.format(dropout_list[idx]))
    print('{:.3f} {:.3f}'.format(accs[idx, 0], accs[idx, 1]))
    print('{:.3f} {:.3f}'.format(stds[idx, 0], stds[idx, 1]))



if __name__=='__main__':
    dataset='OCTAD'
    n_task=2
    # df = read_data(dataset)
    # print(df.head())
    # split_data(df, dataset, n_task=n_task)

    cv_result(dataset, n_task)