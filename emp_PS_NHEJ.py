#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:27:27 2024

@author: corneliasheeran
"""
import sympy as sym
import numpy as np
from sympy import MatrixSymbol, Inverse
import sys
import math
from sympy.matrices import Matrix
from sympy import symbols, latex
from sympy import Matrix, matrix_multiply_elementwise
from sympy import *
import seaborn as sns
import matplotlib.pyplot as plt
import tqdm as tqdm
import itertools
import pandas as pd
from sympy.utilities.lambdify import implemented_function
from sympy import lambdify
from scipy.optimize import fsolve


hval = 0.3
sval = 0.95
eplist = [0.5, 0.9]
Cclist = [100, 500, 1000, 2000, 5000, 1e4]
nulist = np.linspace(0, 0.1, 5)

for epsilonval in eplist:
    for Ccval in tqdm.tqdm(Cclist):
        for nuval in nulist:
            
            round_h  = round(hval, 2)
            round_s  = round(sval, 2)
            round_ep = round(epsilonval, 2)
            round_Cc = round(Ccval, 0)
            round_nu = round(nuval, 2)
            
            format_h  = '{:.2f}'.format(round_h)
            format_s  = '{:.2f}'.format(round_s)
            format_ep = '{:.2f}'.format(round_ep)
            format_Cc = '{:.1f}'.format(round_Cc)
            format_nu = '{:.2f}'.format(round_nu)
                
            path_freqt = f"/rds/general/user/csheeran/home/FREQ_NHEJ/freq_time_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            path_freqk = f"/rds/general/user/csheeran/home/FREQ_NHEJ/freq_space_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            
            path_tD    = f"/rds/general/user/csheeran/home/FFT_t_NHEJ/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            path_tW    = f"/rds/general/user/csheeran/home/FFT_t_NHEJ/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            
            path_kD    = f"/rds/general/user/csheeran/home/FFT_k_NHEJ/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            path_kW    = f"/rds/general/user/csheeran/home/FFT_k_NHEJ/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}_nu={format_nu}.csv"
            
            pd_freqt = pd.read_csv(path_freqt)
            pd_freqk = pd.read_csv(path_freqk)
            pd_tD    = pd.read_csv(path_tD)
            pd_tW    = pd.read_csv(path_tW)
            pd_kD    = pd.read_csv(path_kD)
            pd_kW    = pd.read_csv(path_kW)
            
            freqt = pd_freqt.to_numpy()
            freqk = pd_freqk.to_numpy()
            tD    = pd_tD.to_numpy()
            tW    = pd_tW.to_numpy()
            kD    = pd_kD.to_numpy()
            kW    = pd_kW.to_numpy()
            
        
    
            time = plt.figure(1)
            plt.scatter(freqt, tD, label='Sim Drive')
            plt.scatter(freqt, tW, label='Sim Wild')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.title('Logscale Plot with Scatter Plots')
            plt.xlabel('Temporal Frequency')
            plt.ylabel('Power')
            time.savefig(f"/rds/general/user/csheeran/home/PS_plots/NHEJ_PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            
            spacetime = plt.figure(2)
            plt.scatter(freqk, kD, label='Sim Drive')
            plt.scatter(freqk, kW, label='Sim Wild')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.title('Logscale Plot with Scatter Plots')
            plt.xlabel('Spatial Frequency')
            plt.ylabel('Power')
            spacetime.savefig(f"/rds/general/user/csheeran/home/NHEJ_PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')

