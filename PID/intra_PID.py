import os
import glob
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import gauss_pid as PID
import itertools
#import seaborn as sns
import math
import sys
import random
from tqdm import tqdm


N = 100 #100
sized = 1 #2
num_pairs = 25 #25

list_iter = np.array(range(0, N))

comb = list(itertools.combinations(list_iter, 2))
new_list = []
total = []

for l in list_iter:
   #copy_list = comb
    empty = []
    for k in comb:
        indexes = list(k)
        if l not in indexes:
            indexes.append(l)
            new_list.append(indexes)
            empty.append(indexes)
    total.append(empty)


total = np.array(total) #same target node, pairs of input, what they are 
print(total.shape)
selection = np.random.choice(total.shape[1], size=num_pairs, replace=False)
total = total[:, selection, :] #select 25 random pairs

new_list = np.reshape(total, (total.shape[0]*total.shape[1], total.shape[2]))
print(new_list.shape)

dat_path = "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/Activity_Recordings/" # "/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_L1_mtIII_1k/Activity_Recordings/";
#"/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/Activity_Recordings/"

os.chdir("/home/cs07/Scripts/FUNCTIONS_scripts")
os.chdir('/home/cs07/Scripts/BCT/2019_03_03_BCT')

n = int(sys.argv[1]) #network chosen

data_list_post = []
data_list_pre  = []

trial_list = glob.glob(dat_path + "SN_" + str(n) + "_*.mat")[0]
trial_list = trial_list[0:sized]


rec_data    = sio.loadmat(i)['recorded_data']
post_data_1 = rec_data[:, 7] # MATLAB index starts from 1 but Python starts from 0
post_data   = post_data_1[1:]
data_list_post.append(post_data)
pre_data_1  = rec_data[:, 6] # MATLAB index starts from 1 but Python starts from 0
pre_data    = pre_data_1[1:]
data_list_pre.append(pre_data)

post_data  = np.array(data_list_post)
post_data  = post_data.flatten()
pre_data   = np.array(data_list_pre)
pre_data   = pre_data.flatten()
num_trials = (int(rec_data[-1, 0]) + 1)*sized

M       = post_data.shape[0]
unit    = round(M/(num_trials))
time    = 50
normed  = int(new_list.shape[0]/N)
pid_all = np.zeros((time, normed, N, 4)) #trials, input pairs, target, PID

data_store_post = np.zeros((unit, time, num_trials))
data_store_pre = np.zeros((unit, time, num_trials))
    
for i in range(num_trials):
    k=0
    for j in range((i*unit), ((i+1)*unit)):
        k = k + 1
        data_store_post[k-1, :, i] = np.array(post_data[j])
        data_store_pre[k-1, :, i]  = np.array(pre_data[j])
        

for t in range(0, time):
              
    for g1 in range(total.shape[0]):
        for g2 in range(total.shape[1]):
            
            pid_all[t, g2, g1, :] = PID.mmi(data_store_pre[total[g1, g2, 0], t, :], data_store_pre[total[g1, g2, 1], t, :], data_store_post[total[g1, g2, 2], t, :])


pid_all.tofile("/home/cs07/PIDintra_se/se_" + str(n) + '_.csv', sep = ',')

