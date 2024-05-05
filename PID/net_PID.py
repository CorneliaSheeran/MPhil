import os
import glob
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from Functions import gauss_pid as PID
import itertools
import sys


dat_path = '/Users/corneliasheeran/Documents/*_recording.mat' #"/imaging/shared/users/ja02/CBUActors/SpatialRNNs/Results/mazeGenI_SE1_sWc_mtIV_1k/Activity_Recordings/*.mat"

#os.chdir("/home/cs07/Scripts/FUNCTIONS_scripts")
#import tse_com
#os.chdir('/home/cs07/Scripts/BCT/2019_03_03_BCT')

#sys.path.append("/Users/corneliasheeran/Documents/MPhil_1/scripts/functions")
#import tse_com
#sys.path.append("Users/corneliasheeran/Documents/MPhil_1/BCT/2019_03_03_BCT")

multi = 1;

#def net_PID(task_id):

folder_info = glob.glob(dat_path)[0:multi]
#starting = round((task_id - 1)*100)
#ending   = round(task_id*100 +  1)
n=-1

for file_info in folder_info:
    n=n+1
    rec_data    = sio.loadmat(file_info)['recorded_data']
    post_data_1 = rec_data[:, 7] # MATLAB index starts from 1 but Python starts from 0
    post_data   = post_data_1[1:]

    num_trials  = int(rec_data[-1, 0][0][0])

    M      = post_data.shape[0]
    width  = round(M/(num_trials+1))
    length = len(post_data[4][0])
    
    list_iter = np.array(range(0, 5)) #width
    
    comb = list(itertools.combinations(list_iter, 2))
    new_list = []
    
    for l in list_iter:
        copy_list = comb
        for k in comb:
            indexes = list(k)
            indexes.append(l)
            new_list.append(indexes)
            
    pid_all     = np.zeros((num_trials+1, len(new_list), 4))
    
    for i in range(num_trials+1):
        data_store = np.zeros((width, length))
        
        k = 0
        print(i)
        for j in range((i*width), ((i+1)*width)):
            k = k + 1
            data_store[k-1, :] = np.array(post_data[j][0, :])
            
        g = -1
        for nm in new_list:
            g = g + 1
            pid_all[i, g, :] = PID.mmi(data_store[nm[0], :], data_store[nm[1], :], data_store[nm[2], :])
            
    titled = 'Users/corneliasheeran/Documents/MPhil_1/PID_se_1' #"/home/cs07/PID_se_1"
    plt.plot(np.mean(pid_all, 1)[500:, 0])
    plt.title('PID Synergy Against Trial')
    plt.xlabel('Trial')
    plt.ylabel('Synergy')
    plt.grid(True)
    
    plt.savefig(titled + ".pdf")
