
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 11:50:25 2023

@author: corneliasheeran
"""

import os
import glob
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import gauss_pid as PID
import itertools
import seaborn as sns
import math
import sys
import math
import tqdm as tqdm

dat_path = "/home/cs07/PIDdist_se/"

network = np.arange(0, 1000, 9)


data_list = []
dist_all  = []

file_list = glob.glob(dat_path + "se*.csv")

dist_list = glob.glob(dat_path + "dist*.csv")

for i in tqdm(range(0, network.size)):
    n = network[i]
    file_item = glob.glob(dat_path + f"se_{n}_.csv")
    dist_item = glob.glob(dat_path + f"dist_{n}_.csv")
    
    distance = np.genfromtxt(dist_item, delimiter=',')
    distance = np.reshape(distance, (2500, 3))
    dist_all.append(distance[0:25])
    print(distance[0:5])
    
    data = np.genfromtxt(file_item, delimiter=',')
    data = np.reshape(data, (-1, 2500, 4))
    data_list.append(data[:, 0:25, :])

    

data_arr = np.array(data_list)
data_arr = data_arr[:, :, :, 0]
print(data_arr.shape)
dist_arr = np.array(dist_all)
print(dist_arr.shape)


data_av  = np.nanmean(data_arr, axis=1) #average over trials
print(data_av.shape)

data_av = data_av.flatten()
dist_data = dist_arr[:, :, 0].flatten()


print(data_av.shape)
print(dist_data.shape)

plt.plot(dist_data, data_av)
plt.xlabel('Distance')
plt.ylabel('Synergy')
plt.savefig('/home/cs07/PIDdist_se.pdf')



