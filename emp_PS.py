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
Cclist = [100, 500, 1000, 2000, 5000, 1e4]

for epsilonval in np.linspace(0, 1, 5):
    for Ccval in tqdm.tqdm(Cclist):
            
        round_h  = round(hval, 2)
        round_s  = round(sval, 2)
        round_ep = round(epsilonval, 2)
        round_Cc = round(Ccval, 0)
        
        format_h  = '{:.2f}'.format(round_h)
        format_s  = '{:.2f}'.format(round_s)
        format_ep = '{:.2f}'.format(round_ep)
        format_Cc = '{:.1f}'.format(round_Cc)
            
        path_freqt = f"/rds/general/user/csheeran/home/FREQ/freq_time_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        path_freqk = f"/rds/general/user/csheeran/home/FREQ/freq_space_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        
        path_tD    = f"/rds/general/user/csheeran/home/FFT_t/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        path_tW    = f"/rds/general/user/csheeran/home/FFT_t/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        
        path_kD    = f"/rds/general/user/csheeran/home/FFT_k/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        path_kW    = f"/rds/general/user/csheeran/home/FFT_k/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
        
        #pd_freqt = pd.read_csv(path_freqt)
        #pd_freqk = pd.read_csv(path_freqk)
        #pd_tD    = pd.read_csv(path_tD)
        #pd_tW    = pd.read_csv(path_tW)
        #pd_kD    = pd.read_csv(path_kD)
        #pd_kW    = pd.read_csv(path_kW)
        
        freqt = np.loadtxt(path_freqt, delimiter = ',') #list(pd_freqt.columns)
        freqk = np.loadtxt(path_freqk, delimiter = ',') #list(pd_freqk.columns) #.to_numpy()
        tD    = np.loadtxt(path_tD, delimiter = ',') #list(pd_tD.columns) #.to_numpy()
        tW    = np.loadtxt(path_tW, delimiter = ',') #list(pd_tW.columns) #.to_numpy()
        kD    = np.loadtxt(path_kD, delimiter = ',') #list(pd_kD.columns) #.to_numpy()
        kW    = np.loadtxt(path_kW, delimiter = ',') #list(pd_kW.columns) #.to_numpy()
        #freqt = [float(i) for i in freqt1]
        #freqk = [float(i) for i in freqk]
        #tD = [float(i) for i in tD]
        #tW = [float(i) for i in tW]
        #kD = [float(i) for i in kD]
        #kW = [float(i) for i in kW]

        print(freqt)
       	print("break")
        print(freqk)
        print("break")
        print(tD)
        print("break")
        print(tW)
        print("break")
        print(kD)
        print("break")
        print(kW)
        try:
            time = plt.figure(1)
            plt.scatter(freqt, tD, label='Sim Drive')
            plt.scatter(freqt, tW, label='Sim Wild')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.title('Logscale Temporal Power Sepctra')
            plt.xlabel('Temporal Frequency')
            plt.ylabel('Power')
            time.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            plt.clf()            
        except:
            time = plt.figure(1)
            plt.scatter(freqt, tD, label='Sim Drive')
            plt.scatter(freqt, tW, label='Sim Wild')
            plt.legend()
            plt.title('Temporal Power Sepectra')
            plt.xlabel('Temporal Frequency')
            plt.ylabel('Power')
            time.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            plt.clf()
        try:
            spacetime = plt.figure(2)
            plt.scatter(freqk, kD, label='Sim Drive')
            plt.scatter(freqk, kW, label='Sim Wild')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.title('Logscale Radial Spatial Power Spectra')
            plt.xlabel('Spatial Frequency')
            plt.ylabel('Power')
            spacetime.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_space_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            plt.clf()        
       	except:
            spacetime = plt.figure(2)
            plt.scatter(freqk, kD, label='Sim Drive')
            plt.scatter(freqk, kW, label='Sim Wild')
            plt.legend()
            plt.title('Radial Spatial Power Spectra')
            plt.xlabel('Spatial Frequency')
            plt.ylabel('Power')
            spacetime.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_space_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            plt.clf()
        #time = plt.figure(1)
        #plt.scatter(freqt, tD, label='Sim Drive')
        #plt.scatter(freqt, tW, label='Sim Wild')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.legend()
        #plt.title('Logscale Plot with Scatter Plots')
        #plt.xlabel('Temporal Frequency')
        #plt.ylabel('Power')
        #time.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
        
        #spacetime = plt.figure(2)
        #plt.scatter(freqk, kD, label='Sim Drive')
        #plt.scatter(freqk, kW, label='Sim Wild')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.legend()
        #plt.title('Logscale Plot with Scatter Plots')
        #plt.xlabel('Spatial Frequency')
        #plt.ylabel('Power')
        #spacetime.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')

