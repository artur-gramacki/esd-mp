# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:43:53 2025

@author: Artur 
"""

import h5py
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold

def split_data(hdf5_loc):
    '''
    plit the data into train, validation and test sets
    '''
    
    dir = '../'
    
    names = ["eA_w10_c20_64_64", "eB_w10_c20_64_64", "eC_w10_c20_64_64", "eABC_w10_c20_64_64_XY"]
    
    hdf5_dir = dir + 'working/' + hdf5_loc + "/"
    
    for i in range(4):
        with h5py.File(hdf5_dir + names[i] + '.hdf5',  'r') as f:
            X = f["X"][:]
            Y = f["Y"][:]
            
        x, x_test, y, y_test = train_test_split(X, Y, test_size = 0.1, shuffle = True, random_state = 1024)
        folds = list(StratifiedKFold(n_splits = 5, shuffle = True, random_state = 1024).split(x, y))
        n = len(folds)
        
        del X
        del Y
        
        for fold_number in range(n):
            #print(folds[fold_number][0])
            #print(folds[fold_number][1])
            x_train = x[folds[fold_number][0]]
            y_train = y[folds[fold_number][0]]
            x_validate = x[folds[fold_number][1]]
            y_validate = y[folds[fold_number][1]]
            
            with h5py.File(hdf5_dir + names[i] + '_train' + '_fold_' + str(fold_number) + '.hdf5', 'w') as f:
                f.create_dataset('x_train', data = x_train)
                f.create_dataset('y_train', data = y_train)
            print(names[i] + '_train' + '_fold_' + str(fold_number), flush = True)
                
            with h5py.File(hdf5_dir + names[i] + '_validate' + '_fold_' + str(fold_number) + '.hdf5', 'w') as f:
                f.create_dataset('x_validate', data = x_validate)
                f.create_dataset('y_validate', data = y_validate)
            print(names[i] + '_validate' + '_fold_' + str(fold_number), flush = True)
        
            with h5py.File(hdf5_dir + names[i] + '_test' + '_fold_' + str(fold_number) + '.hdf5', 'w') as f:
                f.create_dataset('x_test', data = x_test)
                f.create_dataset('y_test', data = y_test)    
            print(names[i] + '_test' + '_fold_' + str(fold_number), flush = True)

#############################################################################################
# Run the generators
#############################################################################################  
split_data(hdf5_loc = "hdf5_files")
split_data(hdf5_loc = "hdf5_files_STFT")
