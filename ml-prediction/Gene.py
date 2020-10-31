from __future__ import print_function
import os
import os.path
import numpy as np
import sys
if sys.version_info[0] == 2:
    import cPickle as pickle
else:
    import pickle

import torch.utils.data as data
import random
import pickle as pkl


class Gene(data.Dataset):
    # define train(test)_data, train(test)_label
    # define __getitem__ to return typical data point you want.
    def __init__(self, root, dataset, fold='0', train=True):
        self.train = train  # training set or test set
        self.fold = fold

        fn = root + dataset + '_task_2/fold_' + str(fold) + '.pkl'
        with open(fn, 'rb') as f:
            d = pkl.load(f)
            x_train, y_train_g, y_train_t, x_test, y_test_g, y_test_t = d[0], d[1], d[2], d[3], d[4], d[5]

        # now load the picked numpy arrays
        self.train_data, self.train_labels1, self.train_labels2 = x_train, y_train_g, y_train_t
        self.test_data, self.test_labels1, self.test_labels2 = x_test, y_test_g, y_test_t

        # up or down sampling of train data (hold for now)

    def __getitem__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (feature, target) where target is index of the target class.
        """
        if self.train:
            ft, target1, target2 = self.train_data[index], self.train_labels1[index], self.train_labels2[index]
        else:
            ft, target1, target2 = self.test_data[index], self.test_labels1[index], self.test_labels2[index]

        return ft, target1, target2, index

    def __len__(self):
        if self.train:
            return len(self.train_data)
        else:
            return len(self.test_data)


class Gene3(data.Dataset):
    # define train(test)_data, train(test)_label
    # define __getitem__ to return typical data point you want.
    def __init__(self, root, dataset, fold='0', train=True):
        self.train = train  # training set or test set
        self.fold = fold

        fn = root + dataset + '_task_3/fold_' + str(fold) + '.pkl'
        with open(fn, 'rb') as f:
            d = pkl.load(f)
            x_train, y_train_g, y_train_t, y_train_a, x_test, y_test_g, y_test_t, y_test_a = d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]

        # now load the picked numpy arrays
        self.train_data, self.train_labels1, self.train_labels2, self.train_labels3 = x_train, y_train_g, y_train_t, y_train_a
        self.test_data, self.test_labels1, self.test_labels2, self.test_labels3 = x_test, y_test_g, y_test_t, y_test_a

        # up or down sampling of train data (hold for now)

    def __getitem__(self, index):
        """
        Args:
            index (int): Index

        Returns:
            tuple: (feature, target) where target is index of the target class.
        """
        if self.train:
            ft, target1, target2, target3 = self.train_data[index], self.train_labels1[index], self.train_labels2[index], self.train_labels3[index]
        else:
            ft, target1, target2, target3 = self.test_data[index], self.test_labels1[index], self.test_labels2[index], self.test_labels3[index]

        return ft, target1, target2, target3, index

    def __len__(self):
        if self.train:
            return len(self.train_data)
        else:
            return len(self.test_data)









