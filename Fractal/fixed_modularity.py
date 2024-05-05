#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 09:37:45 2024

@author: corneliasheeran
"""

import math
import numpy as np
import networkx as nx
import random
import scipy.linalg
import scipy.sparse
from scipy.linalg import eigh
import scipy.linalg
import scipy.sparse

# see https://tstojan.github.io/pub/GraphsPrescribedModularity_ComCom.pdf for more information

def D_cal(L_in, C):
    
    L_out = np.sum(C, 0)  #outside links
    D     = 2*L_in + L_out #inside - outside
    
    return D


def delta_m1_cal(L, C, L_in, c, i, j):
    D = D_cal(L_in, C)    
    return (2*L + D[j] - D[i] - 1)/(2*(L**2))



def delta_m2_cal(L, C, L_in, c, i, j):
    D = D_cal(L_in, C)
    return (D[i] - D[j] - 2)/(L**2)


def mincal_delta_m1(L, C, L_in, c):
    
    delta_m1_matrix = np.zeros((c, c))
    
    for i in range(0, c):
        for j in range(0, c):
            
            delta_m1_matrix[i, j] = delta_m1_cal(L, C, L_in, c, i, j)
            
    n, m = np.unravel_index(np.argmin(delta_m1_matrix, axis=None), delta_m1_matrix.shape)
    
    return n, m


def mincal_delta_m2(L, C, L_in, c):
    
    delta_m2_matrix = np.zeros((c, c))
    
    for i in range(0, c):
        for j in range(0, c):
            
            delta_m2_matrix[i, j] = delta_m2_cal(L, C, L_in, c, i, j)
            
    n, m = np.unravel_index(np.argmin(delta_m2_matrix, axis=None), delta_m2_matrix.shape)
    
    return n, m


def Initialise(L, c):
    
    C    = np.zeros((c, c))
    L_in = np.zeros((c,))
    
    r       = L%c 
    k       = math.floor(L/c)
    m_max   = 1 - 1/c - (c-1)/L
    L_in[0] = k
    i       = 2
    
    if r==0:
        while i <= c:
            #create link betwen communities
            C[i-1, i-2] += 1 
            C[i-2, i-1] += 1
            
            L_in[i-1] = k-1
            
            i += 1

        m_max = m_max - 1/(2*(L**2))
        
    elif r <= math.floor(c/2):
        while i <= c-r:
            C[i-1, i-2] += 1 
            C[i-2, i-1] += 1
            
            L_in[i-1] = k-1
        
            if i <= r:
                C[i-1, c-i] += 1 
                C[c-i, i-1] += 1
                
                L_in[c-i] = k
        
            i += 1
        
        L_in[i-1] = k
        m_max = m_max - r*(c-2*r)/(2*c*(L**2))
    
    else:
        while i <= r:
            C[i-1, i-2] += 1 
            C[i-2, i-1] += 1
            
            L_in[i-1] = k
            
            if i <= c-r:
                C[i-1, c-i] += 1 
                C[c-i, i-1] += 1
                
                L_in[i-1] = L_in[i-1] - 1
                L_in[c-i] = k
                
            i += 1
            
        m_max = m_max - (c-r)*(2*r - c)/(2*c*(L**2))
                
    
    return C, m_max, L_in


#state 1 = m is above m_tar
#state 2 = m is below m_tar 

def Replace_Internal_External(C, L, c, m, m_cur, delta_mcur, state, L_in):
    
    i, j = mincal_delta_m1(L, C, L_in, c)
    
    if m_cur > m: #in state 1
    
        delta_m1 = delta_m1_cal(L, C, L_in , c, i, j)
    
        if state == 2 and delta_m1 >= delta_mcur:
            return False, C, L_in, m_cur, delta_mcur, state
        
        if L_in[j-1] == 0:
            return False, C, L_in, m_cur, delta_mcur, state ##
        
        C[i-1, j-1] += 1 
        C[j-1, i-1] += 1
        
        delta_mcur = delta_m1_cal(L, C, L_in , c, i, j)
        m_cur = m_cur - delta_mcur
        
        L_in[j-1] = L_in[j-1] - 1
        state = 1
    
    else: #in state 2
    
        delta_m1 = delta_m1_cal(L, C, L_in , c, i, j)
    
        if state == 1 and delta_m1 >= delta_mcur:
            return False, C, L_in, m_cur, delta_mcur, state  
    
        delta_mcur = delta_m1_cal(L, C, L_in , c, i, j)
        m_cur = m_cur - delta_mcur
        
        if C[i-1, j-1] == 0:
            return False, C, L_in, m_cur, delta_mcur, state ##
        
        if C[i-1, j-1] != 0:
            C[i-1, j-1] = C[i-1, j-1] - 1
            C[j-1, i-1] = C[j-1, i-1] - 1
        else:
            return True, C, L_in, m_cur, delta_mcur, state
        
        L_in[j-1] = L_in[j-1] + 1
        state = 2
    
    return True, C, L_in, m_cur, delta_mcur, state




def Shift_Internal(C, L, c, m, m_cur, delta_mcur, state, L_in):
    
    i, j = mincal_delta_m2(L, C, L_in, c)
    
    if m_cur > m: #in state 1
    
        delta_m2 = delta_m2_cal(L, C, L_in , c, i, j)
    
        if state == 2 and delta_m2 >= delta_mcur:
            return False, C, L_in, m_cur, delta_mcur, state
        
        delta_mcur = delta_m2_cal(L, C, L_in , c, i, j)
        m_cur = m_cur - delta_mcur
        
        L_in[i-1] = L_in[i-1] + 1
        L_in[j-1] = L_in[j-1] - 1
        state = 1
    
    else: #in state 2
    
        delta_m1 = delta_m1_cal(L, C, L_in , c, i, j)
    
        if state == 1 and delta_m1 >= delta_mcur:
            return False, C, L_in, m_cur, delta_mcur, state  
    
        delta_mcur = delta_m2_cal(L, C, L_in , c, i, j)
        m_cur = m_cur - delta_mcur
        
        L_in[i-1] = L_in[i-1] - 1
        L_in[j-1] = L_in[j-1] + 1
        state = 2
    
    return True, C, L_in, m_cur, delta_mcur, state



def rand_trans(p):
    # Generate a random number between 0 and 1
    random_number = random.uniform(0, 1)

    if random_number < p:
        return 1  
    else:
        return 2
    


def TMGG(L, c, m, p, algVariant, epsilon):
    
    C, m_max, L_in = Initialise(L, c)
    m_cur = m_max
    
    if m_cur - epsilon > m:
        print(f'No Graph Can be Generated, for L={L}, c={c} and m={m}')
        C    = float('nan')
        L_in = float('nan')
        return C, L_in
    else:
    
        delta_mcur = +float('inf')
        state = 0
        
        if algVariant == 'StartReplacing':
            
            approachM = True #m_cur < m and
                    
            while abs(m_cur - m) > epsilon and approachM == True: #increases
                if m_cur > 5:
                    print('HELP')
                    break
                else:
                    approachM, C, L_in, m_cur, delta_mcur, state = Replace_Internal_External(C, L, c, m, m_cur, delta_mcur, state, L_in)
                    #print(m_cur)
            approachM = True
                
            while abs(m_cur - m) > epsilon and approachM == True: #decreases
                if m_cur > 5:
                    print('HELP')
                    break
                else:
                    approachM, C, L_in, m_cur, delta_mcur, state = Shift_Internal(C, L, c, m, m_cur, delta_mcur, state, L_in)
                    #print(m_cur)
        return C, L_in    
            
    # elif algVariant == 'StartShifting':
    #     approachM = True
        
    #     while abs(m_cur - m) > epsilon and approachM == True:
    #         approachM = Shift_Internal(C, L, c, m, m_cur, delta_mcur, state, L_in)
             
    #     approachM = True
            
    #     while abs(m_cur - m) > epsilon and approachM == True:
    #         approachM = Replace_Internal_External(C, L, c, m, m_cur, delta_mcur, state, L_in)
            
        
    # elif algVariant == 'Random':
    #     approachM = True
        
    #     while abs(m_cur - m) > epsilon and approachM == True:
            
    #         choice = rand_trans(p)
            
    #         if choice == 1:
    #             approachM = Replace_Internal_External(C, L, c, m, m_cur, delta_mcur, state, L_in)
    #         else:
    #             approachM = ShiftInternal(C, L, c, m, m_cur, delta_mcur, state, L_in)

    # else:
    #     print('Error no algorithm given')
    
    # if approachM == False:
    #     print('There is no graph with modularity in [m-epsilon, m+epsilon]. SECOND')
    # print('m_cur:')
    # print(m_cur)  
    
        #return C, L_in

def fractal_level1(startC, startL, com, mod, epsilon, A_start, totnode):
    
    L_level = np.zeros((com, com)) #old com, new com
    C_level = np.zeros((com, com, com)) #old com, new matrix
    A_out   = np.zeros(np.shape(A_start))
    
    num_total = round(totnode/com)
    num_nodes = round(totnode/com**2)

    if np.isnan(startL).all():
        print('ISNAN')
        return float('nan'), float('nan'), float('nan')
    else:
        for i in range (0, com):
            newl = startL[i]
            C1, L1 = TMGG(newl, com, mod, 0.5,'StartReplacing', epsilon)
            
            L_level[i, :]    = L1
            C_level[i, :, :] = C1
            
            A1 = Random_Graph(L1, num_nodes, com, C1)
            A_out[i*num_total:(i+1)*num_total, i*num_total:(i+1)*num_total] = A1
        
        return C_level, L_level, A_out


def fractal_level2(startC, startL, com, mod, epsilon, A_start, totnode):
    
    L_level = np.zeros((com, com, com)) #oldold com, old com, new com
    C_level = np.zeros((com, com, com, com)) #oldold com, old com, new matrix
    
    num_total = round(totnode/com)
    num2 = round(totnode/com**2)
    num_nodes = round(totnode/com**3)
    A_out   = np.zeros(np.shape(A_start))
    
    if np.isnan(startL).all():
        print('ISNAN')
        return float('nan'), float('nan'), float('nan')
        
    else:
    
        for i in range(0, com):
            for j in range(0, com):
                newl = startL[i, j]
                C1, L1 = TMGG(newl, com, mod, 0.5,'StartReplacing', epsilon)
                A1 = Random_Graph(L1, num_nodes, com, C1)
                
                #A_full[i, j, :, :] = A1
                L_level[i, j, :] = L1
                C_level[i, j, :, :] = C1
                
                A1 = Random_Graph(L1, num_nodes, com, C1)
                A_out[i*num_total + j*num2:i*num_total + (j+1)*num2, i*num_total + j*num2:i*num_total + (j+1)*num2] = A1
        
        
        
        return C_level, L_level, A_out



def Random_Graph(L_in, num_nodes, c, C):
    
    if np.isnan(L_in).all():
        A_full = float('nan')
        return A_full
    else:
    
        A_full = np.zeros((num_nodes*c, num_nodes*c))
        
        for i, n in enumerate(L_in):
            A_n = nx.adjacency_matrix(nx.gnm_random_graph(num_nodes, n, seed=int(42+n), directed=False)).todense()
            A_full[i * num_nodes:(i + 1) * num_nodes, i * num_nodes:(i + 1) * num_nodes] = A_n
        
        
        for i in range(0, c):
            for j in range(0, c):
                if i != j:
                    connect = C[i, j] #number of inter-edges
                    
                    #Ai = A_full[i * num_nodes:(i + 1) * num_nodes, i * num_nodes:(i + 1) * num_nodes]
                    #Aj = A_full[j * num_nodes:(j + 1) * num_nodes, j * num_nodes:(j + 1) * num_nodes]
                    
                    for n in range(0, int(connect)):
                        random_i = random.randint(0, num_nodes-1)
                        random_j = random.randint(0, num_nodes-1)
                        
                        A_full[i * num_nodes + random_i, j * num_nodes + random_j] = 1
                        A_full[j * num_nodes + random_j, i * num_nodes + random_i] = 1
        return A_full

def equationroots(a, b, c): 
 
    # calculating discriminant using formula
    dis = b * b - 4 * a * c 
    sqrt_val = math.sqrt(abs(dis)) 
     
    val1 = (-b + sqrt_val)/(2 * a)
    val2 = (-b - sqrt_val)/(2 * a)
     
    return [val1, val2]

def place_on_diagonal(arrays):
    # Determine the size of the larger array
    total_size = sum(array.shape[0] for array in arrays)
    larger_array = np.zeros((total_size, total_size))

    # Fill the diagonal with the input arrays
    start_row = 0
    for array in arrays:
        size = array.shape[0]
        larger_array[start_row:start_row+size, start_row:start_row+size] = array
        start_row += size

    return larger_array
        


def Complete_Graph(L_in, c, C):
    
    A_list = []
    node_list = []
    
    for i, n in enumerate(L_in):
        
        val_list = equationroots(1, -1, -2*n) 
        nodes    = int(max(val_list))
        
        node_list.append(nodes)
        
        A_n = nx.adjacency_matrix(nx.complete_graph(nodes)).todense()
        A_list.append(A_n)

    A_full  = place_on_diagonal(A_list)
    
    for i in range(0, c):
        for j in range(0, c):
            if i != j:
                connect = C[i, j] #number of inter-edges
                
                sizei = node_list[i]
                sizej = node_list[j]
                
                starti = sum(node_list[:i]) - 1
                startj = sum(node_list[:j]) - 1

                for n in range(0, int(connect)):
                    random_i = random.randint(0, sizei-1)
                    random_j = random.randint(0, sizej-1)
                    
                    A_full[starti + random_i, startj + random_j] = 1
                    A_full[startj + random_j, starti + random_i] = 1
    return A_full



#regular (completely connected) - (complete_graph)
#random (random_graph)
#small world (NWS_Graph)

#hierarchical modular/fractal

#newman_watts_strogatz_graph(n, k, p, seed=None)
#n = int The number of nodes.
#k = int Each node is joined with its k nearest neighbors in a ring topology.
#p = float The probability of adding a new edge for each edge.
#seed = integer, random_state, or None (default)
#edges per community, pick k, k=4, edges = 4*n, nodes=edges/k




def NWS_Graph(L_in, c, C):
    
    A_list    = []
    node_list = []
    
    for i, n in enumerate(L_in):
        

        if n%3 == 0:
            k=3
            nodes = n/3
        elif n%4 == 0:
            k=4
            nodes = n/4
        elif n%5 == 0:
            k=5
            nodes = n/5
        elif n%6 == 0:
            k=6
            nodes = n/6
        elif n%2 == 0:
            k=2
            nodes = n/2
        elif n%1 == 0:
            k=1
            nodes = n/1
            
        node_list.append(int(nodes))
        
        A_n = nx.adjacency_matrix(nx.newman_watts_strogatz_graph(int(nodes), k, 0.25, seed=42)).todense()
        A_list.append(A_n)
    A_full  = place_on_diagonal(A_list)
    #print(node_list)
    for i in range(0, c):
        for j in range(0, c):
            if i != j:
                connect = C[i, j] #number of inter-edges
                
                sizei = node_list[i]
                sizej = node_list[j]
                
                starti = sum(node_list[:i]) - 1
                startj = sum(node_list[:j]) - 1

                for n in range(0, int(connect)):
                    random_i = random.randint(0, sizei-1)
                    random_j = random.randint(0, sizej-1)
                    
                    A_full[starti + random_i, startj + random_j] = 1
                    A_full[startj + random_j, starti + random_i] = 1
    return A_full

#grid_2d_graph(m, n, periodic=False, create_using=None)
#dimensions, m, n int or iterable container of nodes

#A grid graph G_(n,n) has n^2 nodes and 2n^2-2n edges




# def Grid_Graph(L_in, c, C):
    
#     A_list    = []
#     node_list = []
    
#     for i, n in enumerate(L_in):
        
#         val_list = equationroots(1, -1, -n/2)
#         nodes    = int(max(val_list))
            
#         node_list.append(int(nodes))
        
#         A_n = nx.adjacency_matrix(nx.grid_graph(dim=(int(nodes), int(nodes)))).todense()
#         print(A_n.shape)
#         A_list.append(A_n)
        
#     print(node_list)
#     A_full  = place_on_diagonal(A_list)
#     print(A_full.shape)
    
#     for i in range(0, c):
#         for j in range(0, c):
#             if i != j:
#                 connect = C[i, j] #number of inter-edges
                
#                 sizei = node_list[i]
#                 sizej = node_list[j]
                
#                 starti = sum(node_list[:i]) - 1
#                 startj = sum(node_list[:j]) - 1

#                 for n in range(0, int(connect)):
#                     random_i = random.randint(0, sizei-1)
#                     random_j = random.randint(0, sizej-1)
                    
#                     A_full[starti + random_i, startj + random_j] = 1
#                     A_full[startj + random_j, starti + random_i] = 1
#     return A_full



def S(p):

    s = - np.nansum(p * np.log(p))

    return s

def KS_ent(A):
    B = np.copy(A)
    B[B > 0] = 1
    lambda_, _ = eigh(B)
    lambda_max = np.max(np.abs(lambda_))
    H = np.log(lambda_max)
    return H

def Shan_Ent(A_full):
    N = A_full.shape[0]
    D = np.zeros((N, N))
    
    d = np.sum(A_full, axis=1)
    for i in range(N):
        D[i, i] = d[i]
    
    D_inv_sqrt = np.linalg.inv(np.sqrt(D))
    Lag = D - A_full

    M     = np.matmul( np.matmul(D_inv_sqrt, Lag), D_inv_sqrt)
    com   = scipy.linalg.expm(-np.dot(M, 1))  

    S_interim = []
    
    for i in range(0, N):
        E = com[i, :] / np.sum(com[i, :])
        S_interim.append(S(E))

    shan_ent = sum(S_interim)/N
    return shan_ent
