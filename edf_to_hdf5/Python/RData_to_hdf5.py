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

#############################################################################################
# Generate tensors with required dimensions.
# 
# We arrange individual drawings into a 2D matrix
#           ch1 ch2 ch3 ... ch18
# chunk1      x   x   x       x 
# chunk2      x   x   x       x
# .
# .
# .
# chunkN      x   x   x       x
#
#############################################################################################
def create_input_data (mask, image_directory):
    
    # mask: 'fs*' lub 'fns*'
    f = glob.glob(image_directory + mask + '.RData', recursive = False)
    s = int(len(f) / 18) # we have 18 channels, so it's 18
 
    result = pyreadr.read_r(f[0])
    result.keys()
    df = result["tf.matrix"]
    size = df.shape
    
    print("Number of files:", len(f))
    print("Number of samples:", s)
    print("Figure resolution:", size[0], "x", size[1])
    
    # we assume that our images are always square, there is no need 
    # for them to be non-square, it doesn't change anything substantively
    if size[0] != size[1]:
        sys.exit("Png files resolutions must be square.")
    
    n = size[0]
    
    data = np.zeros((n * s, n * 18))
    fileNames = []

    # Sample RData file name: 
    # fns_eA_p03_w10c20_seq_0001_256Hz_ch01_64_64_file00001.RData 
    # fs_eA_p01_w10c20_c01_seq_0001_256Hz_ch01_64_64_file00001.RData
    # File names end with file00001, file00002, etc.
    # For example, files ending with file00001 contain all channels of a given chunk.
    # fs - Filtered Seizure
    # fns - Filtered Non Seizure
    # Filtered - signals from EDF files are subjected to standard filtering:
    #     50Hz notch filter (48.5 Hz - 51.5 Hz)
    #     Low pass IIR Butterworth (30 Hz)
    #     High pass IIR Butterwoth(1 Hz)  
    # p03 - patient number 3
    # w10c20 - window size = 10 and number of contiguous chunks = 20
    # seq_0001 - chunk's next number  
    # c1 - chunk number 1
    # 64_64 - t-f map resolution
    for i in range (0, s):    
        ff = 'file' + str(i + 1).rjust(5, '0')
        files = glob.glob(image_directory + mask + ff + '.Rdata', recursive = False)
        for j in range(0, len(files)):
            img = pyreadr.read_r(files[j])
            df = img["tf.matrix"]
            img = np.float32(df)
            
            if np.isnan(img).any() == True:
                print('NOTE. The file contains NaN values:' + files[j])
                
            data[i * n : (i + 1) * n, j * n : (j + 1) * n] = img
            
            fileNames.append(os.path.basename(files[j]))

        if (i + 1) % 10 == 0:
            print(i + 1, "/", s)
    
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
    
    print('\n')
    print("=====================================================")
    print('X:', X.shape)
    print('Y:', Y.shape)
    print("=====================================================")
    print('\n')

    return(X, Y, fileNames)      
  

dir = '../'

names = ["eA_w10_c20_64_64", "eB_w10_c20_64_64", "eC_w10_c20_64_64"]

for i in range(3):
    image_directory = dir + 'working/tf_maps/' + names[i] + "/"

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
    with h5py.File(dir + 'working/hdf5_files/' + names[i] + '.hdf5', 'w') as f:
        f.create_dataset('X', data = X)
        f.create_dataset('Y', data = Y)
        f.create_dataset('Y_encoded', data = Y_encoded)
        f.create_dataset('fileNames', data = fileNames)  
        
   
with h5py.File(dir + 'working/hdf5_files/' + names[0] + '.hdf5',  'r') as f:
    XA = f["X"][:]
    YA = f["Y"][:]
    
with h5py.File(dir + 'working/hdf5_files/' + names[1] + '.hdf5',  'r') as f:
    XB = f["X"][:]
    YB = f["Y"][:]

with h5py.File(dir + 'working/hdf5_files/' + names[2] + '.hdf5',  'r') as f:
    XC = f["X"][:]
    YC = f["Y"][:]
    
X = np.concatenate((XA, XB, XC))
Y = np.concatenate((YA, YB, YC))

print(X.shape)
print(Y.shape)  

del XA
del XB
del XC
del YA
del YB
del YC  

with h5py.File(dir + 'working/hdf5_files/eABC_w10_c20_64_64_XY.hdf5', 'w') as f:
    f.create_dataset('X', data = X)
    f.create_dataset('Y', data = Y)
    
del X
del Y    

print('done')

#############################################################################################
# Read and plot a sample picture
#############################################################################################
with h5py.File(dir + 'working/hdf5_files/eABC_w10_c20_64_64_XY.hdf5', 'r') as f:
    X = f["X"][:]
    Y = f["Y"][:]

plt.imshow(X[0,0,:,:,:])
plt.imshow(X[1,17,:,:,:])


