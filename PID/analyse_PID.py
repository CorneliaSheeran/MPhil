
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

dat_path = "/home/cs07/PID_se/"

#n = int(sys.argv[1]) #network chosen, averging over output nodes (node chosen (averging over networks))

data_list = []

file_list = glob.glob(dat_path + "se*.csv")
print(file_list)

file_list = file_list[0:2]

for i in tqdm.tqdm(file_list):
    data = np.genfromtxt(i, delimiter=',')
    #print(data.shape[0]/(50*100*4))
    data = np.reshape(data, (-1, 25, 100, 4))
    #data_av = np.nanmean(np.nanmean(data[:, :, :, 0], 1), 1)
    data_list.append(data)
    print(data.shape)

#for i in glob.glob(dat_path + "se*.csv"):
#    print('good')
#    data = np.genfromtxt(i, delimiter=',')
#    try:
#        data = np.reshape(data, (-1, 50, 100, 4))
#        data_av = np.nanmean(np.nanmean(data[:, :, :, 0], axis=1), axis=1)
#        data_list.append(data)
#        print(data.shape)
#    except:
#        print('bad')

data_arr = np.array(data_list)
print(data_arr.shape)


data_av  = np.nanmean(data_arr, axis=0)
print(data_av.shape)
#data_av = np.reshape(data_arr, (data.shape[0], data.shape[1], data.shape[2], data.shape[3]))
#print(data_arr.shape)

#sns.heatmap(data_av)
#plt.savefig("/home/cs07/PID_heatmap/se_heatmap.pdf", dpi=300, bbox_inches='tight')

data_av = np.nanmean(data_av, axis=1)
print(data_av.shape)
#data_av	= np.reshape(data_arr, (data.shape[0], data.shape[2], data.shape[3]))
#print(data_arr.shape)


sns.heatmap(data_av[:, :, 0])
plt.savefig("/home/cs07/PID_heatmap/se_heatmap.pdf", dpi=300, bbox_inches='tight')

data_av = np.nanmean(data_av, axis=1)
print(data_av.shape)
#data_av = np.reshape(data_arr[:, 0], (data.shape[0]))
#print(data_av.shape)

#sns.heatmap(data_av)
#plt.savefig("/home/cs07/PID_heatmap/se_" + str(n) + '_' + 'heatmap.pdf', dpi=300, bbox_inches='tight')

#input_data = glob.glob(dat_path + "se_inputav_" + str(n) + "_.csv")
#data = np.genfromtxt(input_data, delimiter=',')

#for i in range(100):
#    plt.plot(data[:, i])
#    plt.xlabel('Trials')
#    plt.ylabel('Synergy')

#print(data_av)
plt.plot(data_av[:, 0])
plt.xlabel('Trials')
plt.ylabel('Synergy')
plt.savefig('/home/cs07/PIDan_se.pdf')
