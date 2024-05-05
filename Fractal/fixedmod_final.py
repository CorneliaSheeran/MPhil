#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:47:05 2024

@author: corneliasheeran
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from fixed_modularity import *
from networkx import *

    
  
def Get_ent(A_full):

    # try:
        
    if np.isnan(A_full).all():
        print('ISNAN')
        return float('nan'), float('nan')
    else:
        A_full_nx = nx.from_numpy_array(A_full)

        try:
            shan_ent = Shan_Ent(A_full)
            Mod = nx.community.modularity(A_full_nx, nx.community.label_propagation_communities(A_full_nx), weight='weight', resolution=1)           
            
            print('DONE')
            return Mod, shan_ent
        
        except:
            print('SINGULAR')
            return float('nan'), float('nan')
        
    return Mod, shan_ent


def Get_Graph(L, c, epsilon, fixed_mod):
    
        Mod = 0
        ks_ent = 0
        shan_ent = 0
        state = False
        
        # while state == False:
        #     #print(f'state={state}')
            
        C, L_in = TMGG(L, c, fixed_mod, 0.3, 'StartReplacing', epsilon)
        print(L_in)
        
        num_nodes = round(100/c)

            
        A_full = Random_Graph(L_in, num_nodes, c, C) #Random_Graph(L_in, num_nodes, c, C)
        
        #print(A_full)
        if np.isnan(A_full).all():
            return float('nan'), float('nan'), float('nan')
            print('ISNAN')
        else:
            A_full_nx = nx.from_numpy_array(A_full)
            
            try:
                Mod_pre = Mod
                shan_ent = Shan_Ent(A_full)
                Mod = nx.community.modularity(A_full_nx, nx.community.label_propagation_communities(A_full_nx), weight='weight', resolution=1)
                ks_ent   = KS_ent(abs(A_full))            
                
                print('DONE')
                state = True
                return Mod, ks_ent, shan_ent

                        
            except:
                print('SINGULAR')
                return float('nan'), float('nan'), float('nan')

        
ending    = 3                           #number of iterations 
epsilon   = 0.05                        #modularity cut off (mod_fix - epsilon <= mod_true <= mod_fix +epsilon)
com       = 3                           #number of communities 
totnode   = 81                          #total number of nodes in graph
num_nodes = round(totnode/com)          #number of nodes in each community
start_L   = 1000                        #number of edges in whole graph
mod_fix   = np.arange(0.1, 2, 0.1)      #fixed modularity (not equal to true modularity)

Dat_fractal = np.zeros((len(mod_fix), ending, 3, 2))

for i, mod in enumerate(mod_fix): #4-20
    print(f'MOD = {mod}')

    for k in range(0, ending):
        
        C, L = TMGG(start_L, com, mod, 0.5, 'StartReplacing', epsilon)
        A0   = Random_Graph(L, num_nodes, com, C)
        
        C1, L1, A1 = fractal_level1(C, L, com, mod, epsilon, A0, totnode)
        C2, L2, A2 = fractal_level2(C1, L1, com, mod , epsilon, A1, totnode)
        
        mod0, shan0 = Get_ent(A0)
        mod1, shan1 = Get_ent(A1)
        mod2, shan2 = Get_ent(A2)
        
        # print('LOOK HERE')
        # print(np.sum(A0-A1))
        # print(np.sum(A0-A2))
        # print('OKAY DONE')
        
        try:
        
            Dat_fractal[i, k, 0, 0] = mod0
            Dat_fractal[i, k, 0, 1] = shan0 #sigma(from_numpy_array(A0), niter=100, nrand=10, seed=None) #shan0
            Dat_fractal[i, k, 1, 0] = mod1
            Dat_fractal[i, k, 1, 1] = shan1 #sigma(from_numpy_array(A1), niter=100, nrand=10, seed=None)
            Dat_fractal[i, k, 2, 0] = mod2
            Dat_fractal[i, k, 2, 1] = shan2 #sigma(from_numpy_array(A2), niter=100, nrand=10, seed=None) #shan2
        except:
            Dat_fractal[i, k, 0, 0] = math.nan
            Dat_fractal[i, k, 0, 1] = math.nan#sigma(from_numpy_array(A0), niter=100, nrand=10, seed=None) #shan0
            Dat_fractal[i, k, 1, 0] = math.nan#mod1
            Dat_fractal[i, k, 1, 1] = math.nan#sigma(from_numpy_array(A1), niter=100, nrand=10, seed=None)
            Dat_fractal[i, k, 2, 0] = math.nan#mod2
            Dat_fractal[i, k, 2, 1] = math.nan#sigma(from_numpy_array(A2), niter=100, nrand=10, seed=None) #shan2
            

for i in range(0, 3):
    newdat = np.concatenate((Dat_fractal[:, 0, i, :], Dat_fractal[:, 1, i, :], Dat_fractal[:, 2, i, :])) #, Dat_fractal[:, 3, i, :])) #, Dat_fractal[:, 4, i, :])) #, Dat_fractal[:, 5, i, :], Dat_fractal[:, 6, i, :], Dat_fractal[:, 7, i, :], Dat_fractal[:, 8, i, :], Dat_fractal[:, 9, i, :]), 0) # ,Dat_fractal[:, 2, i, :], 1) #, Dat_fractal[:, 3, i, :], Dat_fractal[:, 4, i, :]), 1)
    plt.scatter(newdat[:, 0], newdat[:, 1], label=f'L = {i}')
    
plt.xlabel('Modularity')
plt.ylabel('Shannon Entropy')
plt.title('Random Fractal Graph of Given Level', weight='bold')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
plt.grid(True)
plt.legend()
plt.show()
plt.clf()  



c_fix     = np.arange(3, 20, 1) #number of communities

Data = np.zeros((len(c_fix), len(mod_fix), ending, 2))

for i, fm in enumerate(c_fix): #4-20
    print(f'COM = {fm}')
    for j, moded in enumerate(mod_fix):
        print(f'MOD = {moded}')

        for k in range(0, ending):

            L = 1000
            
            mod_r, ks_ent_r, shan_ent_r       = Get_Graph(L, fm, epsilon, moded) 

            Data[i, j, k, 0] = mod_r
            Data[i, j, k, 1] = shan_ent_r
            
cmap = plt.cm.get_cmap('viridis', 18)

#use DATA_randgraph in git as generating more data (called Data) takes forever.

for i in range(0, 18):
    newdat = np.concatenate((DATA_randgraph[i, :, 0, :], DATA_randgraph[i, :, 1, :], DATA_randgraph[i, :, 2, :], DATA_randgraph[i, :, 3, :], DATA_randgraph[i, :, 4, :]), 1)
    fig = plt.scatter(newdat[:, 0], newdat[:, 1], color=cmap(i), label=f'Coms={i+3}')
    
plt.xlabel('Modularity')
plt.grid(True)
plt.ylabel('Shannon Entropy')
plt.title('Random Graph With Differing Number of Communities', weight='bold')
cbar = plt.colorbar() 
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
#cbar = plt.colorbar(ticks=[-1, 0, 1], labels=['3', '12', '20'])
cbar.set_ticks(ticks=[0, 0.5, 1], labels=['3', '12', '20'])

plt.show()
plt.clf()   

    
