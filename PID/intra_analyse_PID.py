#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 10:58:31 2023

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

dat_path = "/home/cs07/PIDintra_se/"

data_list = []

file_list = glob.glob(dat_path + "se*.csv")

file_list = file_list[0:2]

for i in tqdm.tqdm(file_list):
    data = np.genfromtxt(i, delimiter=',')
    data = np.reshape(data, (50, 25, 100, 4))
    data_list.append(data)
    print(data.shape)


data_arr = np.array(data_list)
print(data_arr.shape)


data_av  = np.nanmean(data_arr[:, :, :, :, 0], axis=0) #average over networks
print(data_av.shape)
#data_av = np.reshape(data_arr, (data.shape[0], data.shape[1], data.shape[2], data.shape[3]))
#print(data_arr.shape)

#sns.heatmap(data_av)
#plt.savefig("/home/cs07/PID_heatmap/se_heatmap.pdf", dpi=300, bbox_inches='tight')

data_av = np.nanmean(data_av, axis=1)
print(data_av.shape)

data_av = np.nanmean(data_av, axis=1)
print(data_av.shape)


plt.plot(data_av)
plt.xlabel('Trials')
plt.ylabel('Synergy')
plt.savefig('/home/cs07/PIDintra_se.pdf')
