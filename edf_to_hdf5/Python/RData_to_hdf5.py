# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:43:53 2025

@author: Artur 
"""

import h5py
import numpy as np
import os
import sys
import glob
import pyreadr  
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold

print(os.getcwd())

def create_input_data(mask, image_directory):
    '''
    Generate tensors with required dimensions.
    '''
    
    # mask: 'fs*' lub 'fns*'
    f = glob.glob(image_directory + mask + '.RData', recursive = False)
    s = int(len(f) / 18) # we have 18 channels, so it's 18
 
    result = pyreadr.read_r(f[0])
    result.keys()
    df = result["tf.map.resampled"]
    size = df.shape
    
    print("Number of files:", len(f), flush = True)
    print("Number of samples:", s, flush = True)
    print("Figure resolution:", size[0], "x", size[1], flush = True)
    
    # we assume that our images are always square, there is no need 
    # for them to be non-square, it doesn't change anything substantively
    if size[0] != size[1]:
        sys.exit("Png files resolutions must be square.")
    
    n = size[0]
    
    data = np.zeros((n * s, n * 18))
    fileNames = []
    
    # File names end with 'file00001', 'file00002', etc.
    # For example, files ending with file00001 contain all channels of a given fragment.
    # fs - Filtered Seizure
    # fns - Filtered Non Seizure
    # Filtered - signals from EDF files are subjected to standard filtering:
    #     50Hz notch filter (48.5 Hz - 51.5 Hz)
    #     Low pass IIR Butterworth (30 Hz)
    #     High pass IIR Butterwoth(1 Hz)  
    for i in range (0, s):    
        ff = 'file' + str(i + 1).rjust(5, '0')
        files = glob.glob(image_directory + mask + ff + '.RData', recursive = False)
        for j in range(0, len(files)):
            img = pyreadr.read_r(files[j])
            df = img["tf.map.resampled"]
            img = np.float32(df)
            
            if np.isnan(img).any() == True:
                print('NOTE. The file contains NaN values:' + files[j])
                
            data[i * n : (i + 1) * n, j * n : (j + 1) * n] = img
            
            fileNames.append(os.path.basename(files[j]))

        if (i + 1) % 10 == 0:
            print(i + 1, "/", s, flush = True)
    
    X = np.zeros((s, 18, n, n))
        
    for i in range(0, s):
        for j in range(0, 18):
            X[i,j,:,:] = data[i * n : (i + 1) * n, j * n: (j + 1) * n]

    if mask == 'fs*':
        Y = np.ones(s, dtype = np.byte)
        
    if mask == 'fns*':
        Y = np.zeros(s, dtype = np.byte)       
        
    # For conv3D the last dimension must be 'channels'
    # https://www.tensorflow.org/api_docs/python/tf/keras/layers/Conv3D
    # If data_format="channels_last": 5D tensor with shape: 
    # (batch_size, spatial_dim1, spatial_dim2, spatial_dim3, channels)
    X = np.expand_dims(X, axis = 4)        
    
    print('\n', flush = True)
    print("=====================================================", flush = True)
    print('X:', X.shape, flush = True)
    print('Y:', Y.shape, flush = True)
    print("=====================================================", flush = True)
    print('\n', flush = True)

    return(X, Y, fileNames)      

def create_hdf5_files(hdf5_loc, tf_maps_loc):
    '''
    Generate hdf5 files.
    '''
    
    dir = '../'
    
    names = ["eA_w10_c20_64_64", "eB_w10_c20_64_64", "eC_w10_c20_64_64"]
    
    for i in range(3):
        image_directory = dir + 'working/' + tf_maps_loc + "/" + names[i] + "/"
    
        x_fs, y_fs, fnames_fs = create_input_data(mask = 'fs*', image_directory = image_directory)
        x_fns, y_fns, fnames_fns = create_input_data(mask = 'fns*', image_directory = image_directory)
    
        X = np.concatenate((x_fns, x_fs))
        Y = np.concatenate((y_fns, y_fs))
        fileNames = fnames_fns + fnames_fs
        Y_encoded = np.concatenate((y_fns, y_fs))
    
        # Converting Y to one-hot-encoding
        Y_encoded = np.ones((Y.size, Y.max() + 1), dtype = np.byte)
        Y_encoded[np.arange(Y.size), Y] = 0 
    
    
        # Saving to hdf5
        with h5py.File(dir + 'working/' + hdf5_loc + "/" +  names[i] + '.hdf5', 'w') as f:
            f.create_dataset('X', data = X)
            f.create_dataset('Y', data = Y)
            f.create_dataset('Y_encoded', data = Y_encoded)
            f.create_dataset('fileNames', data = fileNames)  
            
       
    with h5py.File(dir + 'working/' + hdf5_loc + "/" + names[0] + '.hdf5',  'r') as f:
        XA = f["X"][:]
        YA = f["Y"][:]
        
    with h5py.File(dir + 'working/' + hdf5_loc + "/" + names[1] + '.hdf5',  'r') as f:
        XB = f["X"][:]
        YB = f["Y"][:]
    
    with h5py.File(dir + 'working/' + hdf5_loc + "/" + names[2] + '.hdf5',  'r') as f:
        XC = f["X"][:]
        YC = f["Y"][:]
        
    X = np.concatenate((XA, XB, XC))
    Y = np.concatenate((YA, YB, YC))
    
    print(X.shape, flush = True)
    print(Y.shape, flush = True)  
    
    del XA
    del XB
    del XC
    del YA
    del YB
    del YC  
    
    with h5py.File(dir + 'working/' + hdf5_loc + "/" + 'eABC_w10_c20_64_64_XY.hdf5', 'w') as f:
        f.create_dataset('X', data = X)
        f.create_dataset('Y', data = Y)
        
    del X
    del Y    
    
    print('done', flush = True)

#############################################################################################
# Run the generators
#############################################################################################    
# for generating hdf5 files based on Matching Pursuit    
create_hdf5_files(hdf5_loc = "hdf5_files", tf_maps_loc = "tf_maps")

# for generating hdf5 files based on Short Time Fourier Transform
create_hdf5_files(hdf5_loc = "hdf5_files_STFT", tf_maps_loc = "tf_maps_STFT")


#############################################################################################
# Read and plot a sample picture
#############################################################################################
with h5py.File("../" + 'working/' + "hdf5_files/" +  'eABC_w10_c20_64_64_XY.hdf5', 'r') as f:
    X_MP = f["X"][:]
    Y_MP = f["Y"][:]

plt.imshow(X_MP[0,0,:,:,:])
plt.imshow(X_MP[195,17,:,:,:]) 


with h5py.File("../" + 'working/' + "hdf5_files_STFT/" +  'eABC_w10_c20_64_64_XY.hdf5', 'r') as f:
    X_STFT = f["X"][:]
    Y_STFT = f["Y"][:]

plt.imshow(X_STFT[0,0,:,:,:])
plt.imshow(X_STFT[195,17,:,:,:]) 

