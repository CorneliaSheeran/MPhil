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


Bb, h, s, epsilon, V, Cc, Dm, nu, hN = sym.symbols("Bb h s epsilon V Cc Dm nu hN", real=True)
a, ad, b, bd, c, cd, e, ed = sym.symbols("a ad b bd c cd e ed", real=True)

def sp_partial_inversion(m, *cells):
    ''' Partial inversion algorithm.
    ARGUMENTS
        m <sympy.Matrix> : symbolic matrix
        *cells <tuple>   : 2-tuples with matrix indices to perform partial inversion on.
    RETURNS
        <sympy.Matrix>   : matrix m partially-inverted at indices *cells
    '''
    # Partial inversion algorithm
    M = m.copy()
    for cell in cells:
        i,k = cell
        z = M[i,k]
        newmat = []
        for r in range(m.shape[0]):
            newrow = []
            for c in range(m.shape[1]):
                if i == r:
                    if k == c:  # Pivot cell
                        newrow.append( 1/z )
                    else:       # Pivot row
                        newrow.append( -M[r,c]/z )
                else:
                    if k == c:  # Pivot col
                        newrow.append( M[r,c]/z )
                    else:       # Elsewhere
                        newrow.append( M[r,c] - M[i,c]*M[r,k]/z )
            newmat.append(newrow)
        M =  Matrix(newmat)
    #
    return M

def FullInversion(m):
    "Full matrix inversion is partial inversion at all i==j."
    cells = [(i,i) for i in range(m.shape[0])]
    return sp_partial_inversion(m, *cells)


def calc_fit(I1):
    
    if (I1[0] == ad or I1[1] == ad): #one W allele
        if (I1[0] == bd or I1[1] == bd): #WD
            fit = 1-h*s
        elif (I1[0] == cd or I1[1] == cd): #WN
            fit = 1-hN*s
        else:
            fit = 1 #WW
    else:
        fit = 1-s #NN, DN, DD
        
    return fit


def calc_rate(I1, I2, Out):
    
    R = []

    av_fit = calc_fit(I1)
    
    for O1 in Out:
        
        if O1[0] == cd or O1[1] == cd:
            rate = sym.simplify(av_fit/4)
            R.append(rate)
        else:
            
            if I1[0] == cd or I1[1] == cd: #1 N allele
                multi1 = 1
            elif I1[0] == bd or I1[1] == bd: #if one input allele is drive
    
                if I1[0] == bd and I1[1] == bd: #if both input alleles are drive
                    multi1 = 1
                      
                else:
                    if O1[0] == bd:
                        multi1 = (1+epsilon)
                    else: 
                        multi1 = (1-epsilon)
            else:
                multi1 = 1
                
                
            if I2[0] == cd or I2[1] == cd: #1 N allele
                multi2 = 1
                
            elif I2[0] == bd or I2[1] == bd:
    
                if I2[0] == bd and I2[1] == bd:
                    multi2 = 1
                      
                else:
                    if O1[1] == bd:
                        multi2 = (1+epsilon)
                    else:
                        multi2 = (1-epsilon)
            else:
                multi2 = 1
                        
            multi = multi1*multi2
            rate  = sym.simplify(av_fit*multi/4)
            R.append(rate)
        
    return R, av_fit


def H_part(I1, I2, Ibar1, Ibar2, Out):
    H = 0
    Rate, av_fit = calc_rate(I1, I2, Out)
    
    for i in range(0, len(Out)):
        O1 = Out[i]
        R  = Rate[i]
        R4 = sym.simplify(R*4)
        
        if '(1 + epsilon)' in str(R4) or '(epsilon + 1)' in str(R4) or 'epsilon + 1' in str(R4) or '1 + epsilon' in str(R4):
            R  = (R/(1+epsilon))*(1+epsilon*(1-nu))
            H += R*(O1[0]*O1[1] - ed)*(I1[0]*I1[1])*(I2[0]*I2[1])*(Ibar1[0]*Ibar1[1])*(Ibar2[0]*Ibar2[1])*e
            
            R1 = av_fit*epsilon*nu/4
            
            if (O1[0] == bd and O1[1] == bd):
                O2 = (bd, cd)
            else:
                O2 = (ad, cd)
                
            H += R1*(O2[0]*O2[1] - ed)*(I1[0]*I1[1])*(I2[0]*I2[1])*(Ibar1[0]*Ibar1[1])*(Ibar2[0]*Ibar2[1])*e
        
        else:
            H += R*(O1[0]*O1[1] - ed)*(I1[0]*I1[1])*(I2[0]*I2[1])*(Ibar1[0]*Ibar1[1])*(Ibar2[0]*Ibar2[1])*e

    return Bb*H/(V**4)


def H_void(geno_crea, geno_anih):
    H_death = 0
    
    for i in range(0, len(geno_crea)):
        X    = geno_crea[i]
        Xbar = geno_anih[i]
        H_death += ed*(Xbar[0]*Xbar[1]) - X[0]*X[1]*Xbar[0]*Xbar[1]
    return H_death/V


crea = [ad, bd, cd]
anih = [a, b, c]

geno_crea = list(itertools.combinations(crea,2)) + list(zip(crea, crea))
geno_anih = list(itertools.combinations(anih,2)) + list(zip(anih, anih))


H_death = H_void(geno_crea, geno_anih)

#sym.pprint(H_death)
Ht  = -H_death*Dm

for i in range(0, len(geno_crea)):
    for j in range(0, len(geno_crea)):
        
        I1      = geno_crea[i]
        #print(I1)
        I2      = geno_crea[j]
        #print(I2)
        Ibar1   = geno_anih[i]
        Ibar2   = geno_anih[j]
        Out     = [(x, y) for x, y in itertools.product(I1, I2)]
        #print(Out)
        Outbar  = [(x, y) for x, y in itertools.product(Ibar1, Ibar2)]
        
        #print(H_part(I1, I2, Ibar1, Ibar2, Out))
        Ht += -H_part(I1, I2, Ibar1, Ibar2, Out)
        
        #print('\n')

#sym.pprint(sym.diff(H_full, bd))

aCH, adCH, bCH, bdCH, cCH, cdCH, eCH, edCH = sym.symbols("aCH adCH bCH bdCH cCH cdCH eCH edCH", real=True)
z, Zm, Z, y, Ym, Y, x, Xm, X, t, Tm, T = sym.symbols("z Zm Z y Ym Y x Xm X t Tm T", real=True)

aCH  = z*Zm
adCH = Z
bCH  = y*Ym
bdCH = Y
cCH  = x*Xm
cdCH = X
eCH  = t*Tm
edCH = T

zd, yd, xd, td = sym.symbols("zd yd xd td", real=True)
W, D, N, Es    = sym.symbols("W D N Es", real=True)
n1, n2, n3, n4 = sym.symbols("n1 n2 n3 n4", real=True)

solz = V*W  + (V**0.5)*n1
soly = V*D  + (V**0.5)*n2
solx = V*N  + (V**0.5)*n3
solt = V*Es + (V**0.5)*n4


solZ  = 1 + zd/(V**0.5) + (zd**2)/(2*V)
solZm = 1 - zd/(V**0.5) + (zd**2)/(2*V)
solY  = 1 + yd/(V**0.5) + (yd**2)/(2*V)
solYm = 1 - yd/(V**0.5) + (yd**2)/(2*V)
solX  = 1 + xd/(V**0.5) + (xd**2)/(2*V)
solXm = 1 - xd/(V**0.5) + (xd**2)/(2*V)
solT  = 1 + td/(V**0.5) + (td**2)/(2*V)
solTm = 1 - td/(V**0.5) + (td**2)/(2*V)

Ht = Ht.subs({a: aCH, ad: adCH, b: bCH, bd: bdCH, c:cCH, cd: cdCH, e: eCH, ed: edCH})
Ht = Ht.subs({Z*z*Zm: z-1, Y*y*Ym: y-1, X*x*Xm: x-1, T*t*Tm: t-1})
Ht = Ht.subs({Z: solZ, Zm: solZm, Y: solY, Ym: solYm, X: solX, Xm: solXm, T:solT, Tm:solTm})
Ht = Ht.subs({z: solz, y:soly, x:solx, t:solt})

Ht_MF = Ht/(V**0.5)
Ht_PS = Ht
Ht_ex = Ht*(V**0.5)

print('done0')

#L_MF  = sym.collect(sym.expand(Ht_MF), V)
L_PS  = sym.collect(sym.expand(Ht_PS), V)
# L_ex  = sym.collect(sym.expand(Ht_ex), V)
#Ht_MF = L_MF.coeff(V, 0)
Ht_PS = L_PS.coeff(V, 0)
# Ht_ex = L_ex.coeff(V, 0)

print('done1')

Wsol = (-1)*((-Bb*D**3*Es*W*epsilon*h*s/2 - Bb*D**3*Es*W*epsilon*s/2 + Bb*D**3*Es*W*epsilon 
+ Bb*D**3*Es*W*h*s/2 + Bb*D**3*Es*W*s/2 - Bb*D**3*Es*W - Bb*D**2*Es*N*W*epsilon*h*s/4 
- Bb*D**2*Es*N*W*epsilon*s/4 + Bb*D**2*Es*N*W*epsilon/2 + Bb*D**2*Es*N*W*h*s/2 
+ Bb*D**2*Es*N*W*hN*s/2 + Bb*D**2*Es*N*W*s - 2*Bb*D**2*Es*N*W + Bb*D**2*Es*W**2*epsilon**2*h*nu*s/2 
- Bb*D**2*Es*W**2*epsilon**2*nu/2 - Bb*D**2*Es*W**2*epsilon*h*s + Bb*D**2*Es*W**2*epsilon 
+ Bb*D**2*Es*W**2*h*s + Bb*D**2*Es*W**2*s - 3.0*Bb*D**2*Es*W**2 + Bb*D*Es*N**2*W*h*s/2 
+ Bb*D*Es*N**2*W*hN*s/2 + Bb*D*Es*N**2*W*s - 2*Bb*D*Es*N**2*W - Bb*D*Es*N*W**2*epsilon*h*s/4 
- Bb*D*Es*N*W**2*epsilon*hN*s/4 + Bb*D*Es*N*W**2*epsilon/2 + Bb*D*Es*N*W**2*h*s 
+ Bb*D*Es*N*W**2*hN*s + Bb*D*Es*N*W**2*s - 4.0*Bb*D*Es*N*W**2 
+ Bb*D*Es*W**3*epsilon**2*nu/(2*epsilon + 2) - Bb*D*Es*W**3*epsilon**2/(2*epsilon + 2) 
- Bb*D*Es*W**3*epsilon*h*s/2 - Bb*D*Es*W**3*epsilon*nu/2 + Bb*D*Es*W**3*epsilon*nu/(2*epsilon + 2) 
+ 3*Bb*D*Es*W**3*epsilon/2 - Bb*D*Es*W**3*epsilon/(epsilon + 1) + 3*Bb*D*Es*W**3*h*s/2 
- 5*Bb*D*Es*W**3/2 - Bb*D*Es*W**3/(2*epsilon + 2) + Bb*Es*N**3*W*hN*s/2 + Bb*Es*N**3*W*s/2 
- Bb*Es*N**3*W + Bb*Es*N**2*W**2*hN*s + Bb*Es*N**2*W**2*s - 3.0*Bb*Es*N**2*W**2 
+ 3*Bb*Es*N*W**3*hN*s/2 - 3.0*Bb*Es*N*W**3 - 2.0*Bb*Es*W**4)/(2*Cc**3) + (D*Dm*W + Dm*N*W + 2*Dm*W**2)/(2*Cc)) 


Dsol = (-1)*((2*Bb*D**4*Es*s - 2*Bb*D**4*Es + 3*Bb*D**3*Es*N*s - 3*Bb*D**3*Es*N - Bb*D**3*Es*W*epsilon*h*nu*s/2 
+ Bb*D**3*Es*W*epsilon*h*s/2 - Bb*D**3*Es*W*epsilon*nu*s/2 + Bb*D**3*Es*W*epsilon*nu 
+ Bb*D**3*Es*W*epsilon*s/2 - Bb*D**3*Es*W*epsilon + 3*Bb*D**3*Es*W*h*s/2 + 3*Bb*D**3*Es*W*s/2 
- 3*Bb*D**3*Es*W + 3*Bb*D**2*Es*N**2*s - 3*Bb*D**2*Es*N**2 - Bb*D**2*Es*N*W*epsilon*h*nu*s/4 
+ Bb*D**2*Es*N*W*epsilon*h*s/4 - Bb*D**2*Es*N*W*epsilon*nu*s/4 + Bb*D**2*Es*N*W*epsilon*nu/2 
+ Bb*D**2*Es*N*W*epsilon*s/4 - Bb*D**2*Es*N*W*epsilon/2 + Bb*D**2*Es*N*W*h*s + Bb*D**2*Es*N*W*hN*s 
+ 2*Bb*D**2*Es*N*W*s - 4*Bb*D**2*Es*N*W - 3*Bb*D**2*Es*W**2*epsilon*h*nu*s/4 + Bb*D**2*Es*W**2*epsilon*h*s 
+ 3*Bb*D**2*Es*W**2*epsilon*nu/4 - Bb*D**2*Es*W**2*epsilon + Bb*D**2*Es*W**2*h*s + Bb*D**2*Es*W**2*s 
- 3.0*Bb*D**2*Es*W**2 + Bb*D*Es*N**3*s - Bb*D*Es*N**3 + Bb*D*Es*N**2*W*h*s/2 + Bb*D*Es*N**2*W*hN*s/2 
+ Bb*D*Es*N**2*W*s - 2*Bb*D*Es*N**2*W - Bb*D*Es*N*W**2*epsilon*h*nu*s/4 + Bb*D*Es*N*W**2*epsilon*h*s/4 
- Bb*D*Es*N*W**2*epsilon*hN*nu*s/4 + Bb*D*Es*N*W**2*epsilon*hN*s/4 + Bb*D*Es*N*W**2*epsilon*nu/2 
- Bb*D*Es*N*W**2*epsilon/2 + Bb*D*Es*N*W**2*h*s/2 + Bb*D*Es*N*W**2*hN*s/2 + Bb*D*Es*N*W**2*s/2 
- 2.0*Bb*D*Es*N*W**2 + Bb*D*Es*W**3*epsilon**2*nu/(2*epsilon + 2) 
- Bb*D*Es*W**3*epsilon**2/(2*epsilon + 2) - Bb*D*Es*W**3*epsilon*h*nu*s/2 
+ Bb*D*Es*W**3*epsilon*h*s/2 + Bb*D*Es*W**3*epsilon*nu/2 + Bb*D*Es*W**3*epsilon*nu/(2*epsilon + 2) 
- Bb*D*Es*W**3*epsilon/2 - Bb*D*Es*W**3*epsilon/(epsilon + 1) + Bb*D*Es*W**3*h*s/2 - Bb*D*Es*W**3/2 
- Bb*D*Es*W**3/(2*epsilon + 2))/(2*Cc**3) + (2*D**2*Dm + D*Dm*N + D*Dm*W)/(2*Cc)) 

Nsol = (-1)*((Bb*D**3*Es*N*s - Bb*D**3*Es*N + Bb*D**3*Es*W*epsilon*h*nu*s/2 + Bb*D**3*Es*W*epsilon*nu*s/2 
- Bb*D**3*Es*W*epsilon*nu + 3*Bb*D**2*Es*N**2*s - 3*Bb*D**2*Es*N**2 + Bb*D**2*Es*N*W*epsilon*h*nu*s/4 
+ Bb*D**2*Es*N*W*epsilon*nu*s/4 - Bb*D**2*Es*N*W*epsilon*nu/2 + Bb*D**2*Es*N*W*h*s/2 + Bb*D**2*Es*N*W*hN*s/2 
+ Bb*D**2*Es*N*W*s - 2*Bb*D**2*Es*N*W + 3*Bb*D**2*Es*W**2*epsilon*h*nu*s/4 - 3*Bb*D**2*Es*W**2*epsilon*nu/4 
+ 3*Bb*D*Es*N**3*s - 3*Bb*D*Es*N**3 + Bb*D*Es*N**2*W*h*s + Bb*D*Es*N**2*W*hN*s + 2*Bb*D*Es*N**2*W*s 
- 4*Bb*D*Es*N**2*W + Bb*D*Es*N*W**2*epsilon*h*nu*s/4 + Bb*D*Es*N*W**2*epsilon*hN*nu*s/4 
- Bb*D*Es*N*W**2*epsilon*nu/2 + Bb*D*Es*N*W**2*h*s/2 + Bb*D*Es*N*W**2*hN*s/2 + Bb*D*Es*N*W**2*s/2 
- 2.0*Bb*D*Es*N*W**2 + Bb*D*Es*W**3*epsilon*h*nu*s/2 - Bb*D*Es*W**3*epsilon*nu + 2*Bb*Es*N**4*s - 2*Bb*Es*N**4 
+ 3*Bb*Es*N**3*W*hN*s/2 + 3*Bb*Es*N**3*W*s/2 - 3*Bb*Es*N**3*W + Bb*Es*N**2*W**2*hN*s + Bb*Es*N**2*W**2*s 
- 3.0*Bb*Es*N**2*W**2 + Bb*Es*N*W**3*hN*s/2 - 1.0*Bb*Es*N*W**3)/(2*Cc**3) + (D*Dm*N + 2*Dm*N**2 + Dm*N*W)/(2*Cc)) 

x = [n1, n2]
y = [zd, yd]
source = [W, D]


A = Matrix.zeros(2, 2)
B = Matrix.zeros(2, 2)

item, kx, Dx, omega, g1, g2 = sym.symbols("item kx Dx omega g1 g2", real=True)

for i, n in enumerate(y):
    for j, m in enumerate(x):
        inter1  = Ht_PS.subs({m*n: item})
        inter2  = sym.collect(inter1, item)
        inter4  = inter2.coeff(item, 1)
        A[i, j] = inter4    
        if i == j:
            A[i, j] += -kx*Dx

#sym.pprint(sym.simplify(A))

for i, n in enumerate(y):
    for j, m in enumerate(y):
        inter1  = Ht_PS.subs({m*n: item})
        inter2  = sym.collect(inter1, item)
        inter4  = inter2.coeff(item, 1)
        B[i, j] = inter4
        if i == j:
            B[i, j] = 2*inter4
            B[i, j] += source[i]*Dx*kx

print('done2')

Wsol = Wsol.subs({Es:(1-(W+D+N)/(2*Cc))/2})
Wsol = sym.simplify(Wsol.subs({N:0, nu:0, Bb:40, Dm:2}))

Dsol = Dsol.subs({Es:(1-(W+D+N)/(2*Cc))/2})
Dsol = sym.simplify(Dsol.subs({N:0, nu:0, Bb:40, Dm:2}))

A = sym.simplify(A.subs({N:0, nu:0, Bb:40, Dm:2}))
B = sym.simplify(B.subs({N:0, nu:0, Bb:40, Dm:2}))

print('done3')
hval = 0.3
sval = 0.95

#epsilonvaled = np.linspace(0, 1, 5)
#index = int(sys.argv[1])

epsilonval = 0.7 #epsilonvaled[index]
Cclist = [100, 500, 1000, 5000, 1e4]
# for hval in tqdm.tqdm(np.linspace(0, 1, 5)):
#     for sval in np.linspace(0, 1, 5):
#for epsilonval in np.linspace(0, 1, 5):

for Ccval in tqdm.tqdm(Cclist):
    Wsol = sym.expand(Wsol.subs({h:hval, s:sval, epsilon:epsilonval, Cc:Ccval})) 
    Dsol = sym.expand(Dsol.subs({h:hval, s:sval, epsilon:epsilonval, Cc:Ccval}))
    
    equations = [Wsol, Dsol]

    # Create lambda functions for numerical evaluation
    lam_fun = lambdify((W, D), equations)
    
    
    # fun = Wsol + Dsol
    # xvec = sym.Array([W, D])
    # lam_fun = lambdify(xvec, fun)
    
    def to_solve(vars):
        W_val, D_val = vars
        return lam_fun(W_val, D_val)

    initial_guess = [Ccval/2, Ccval/2]
    
    # Solve the system of equations
    Ans = fsolve(to_solve, initial_guess)
    out = to_solve(Ans)
    print(Ans)    
    if out[0] > 1 or out[1] > 1 or floor(Ans[0]) < 0 or floor(Ans[1]) < 0:
        print('stop')
        print(Ccval)
    else: 
        if floor(Ans[0]) < 0 or floor(Ans[1]) < 0:
    
            print('perturbation')
            print(Ccval)       
        else:
            
            print('normal')
            
            A = sym.simplify(A.subs({h:hval, s:sval, epsilon:epsilonval, Cc:Ccval, W:Ans[0], D:Ans[1], Es: (1-(Ans[0] + Ans[1])/(2*Cc))/2, Dx:0.05}))
            B = sym.simplify(B.subs({h:hval, s:sval, epsilon:epsilonval, Cc:Ccval, W:Ans[0], D:Ans[1], Es: (1-(Ans[0] + Ans[1])/(2*Cc))/2, Dx:0.05}))
            print('done')
            shiftA   = sym.eye(2)*sym.I*omega
            Ashifted = A + shiftA
            Dinv     = FullInversion(Ashifted)
            
            print('done')
            
            Gamma = sym.Matrix([[g1], [g2]])
            
            newx    = Dinv*Gamma
            first   = newx
            second  = sym.conjugate(newx)
            PS_full = matrix_multiply_elementwise(first, second)
            
            print('done')
            
            PS_W = sym.expand(PS_full[0])
            PS_D = sym.expand(PS_full[1])
            
            PS_W = sym.expand(PS_W.subs({g1*g1: B[0, 0], g1*g2: B[0, 1], g2*g2:B[1, 1]}))
            PS_D = sym.expand(PS_D.subs({g1*g1: B[0, 0], g1*g2: B[0, 1], g2*g2:B[1, 1]}))

            PS_W = PS_W.as_real_imag()[0]
            PS_D = PS_D.as_real_imag()[0]
            print('done fin')
            #print(sym.simplify(PS_W.as_real_imag()[0]))
            #print(sym.simplify(PS_D.as_real_imag()[0]))
            
            #for temporal-radial power spectral
            
            PS_W_TR = np.zeros((10, 10))
            PS_D_TR = np.zeros((10, 10))
            
            for i, ki in enumerate(np.linspace(0.1, 1, 10)):
                for j, oi in enumerate(np.linspace(0.1, 1, 10)):
                    #print(sym.simplify(sym.expand(PS_W.subs({omega:oi, kx:ki, Dx:0.05}))))
                    try:
                        PS_W_TR[i, j] = sym.simplify(sym.expand(PS_W.subs({omega:oi, kx:ki, Dx:0.05})))
                        PS_D_TR[i, j] = sym.simplify(sym.expand(PS_D.subs({omega:oi, kx:ki, Dx:0.05})))
                    except:
                        print('HELP')
                       
            print(PS_W_TR)
            #for temporal power spectra
            ki     = np.linspace(0.1, 1, 10)
            deltak = ki[1] - ki[0]
    
            PS_W_t = 0
            PS_D_t = 0
      
            for val in ki:
                W_t = PS_W.subs({kx: val})
                D_t = PS_D.subs({kx: val})
                
                PS_W_t += W_t*deltak
                PS_D_t += D_t*deltak
     
            PS_W_t = sym.simplify(sym.expand(PS_W_t.as_real_imag()[0]))
            funPS_W_t = sym.lambdify((omega), PS_W_t, 'numpy')
            PS_D_t = sym.simplify(sym.expand(PS_D_t.as_real_imag()[0]))
            funPS_D_t = sym.lambdify((omega), PS_D_t, 'numpy')
            
            print('done')
            round_h  = round(hval, 2)
            round_s  = round(sval, 2)
            round_ep = round(epsilonval, 2)
            round_Cc = round(Ccval, 0)
            if round_Cc == 762 or round_Cc == 765:
                round_Cc = 763
            
            format_h  = '{:.2f}'.format(round_h)
            format_s  = '{:.2f}'.format(round_s)
            format_ep = '{:.2f}'.format(round_ep)
            format_Cc = '{:.1f}'.format(round_Cc)
                
            np.save(f"/rds/general/user/csheeran/home/sim_dat/simDTR_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.npy", PS_D_TR)
            np.save(f"/rds/general/user/csheeran/home/sim_dat/simWTR_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.npy", PS_W_TR)

            path_freqt = f"/rds/general/user/csheeran/home/FREQ/freq_time_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            path_freqk = f"/rds/general/user/csheeran/home/FREQ/freq_space_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            
            path_tD    = f"/rds/general/user/csheeran/home/FFT_t/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            path_tW    = f"/rds/general/user/csheeran/home/FFT_t/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            
            path_kD    = f"/rds/general/user/csheeran/home/FFT_k/drive_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            path_kW    = f"/rds/general/user/csheeran/home/FFT_k/wild_h={format_h}_s={format_s}_ep={format_ep}_Cc={format_Cc}.csv"
            
            freqt = np.genfromtxt(path_freqt, delimiter=',')
            freqk = np.genfromtxt(path_freqk, delimiter=',')
            tD    = np.genfromtxt(path_tD, delimiter=',')
            tW    = np.genfromtxt(path_tW, delimiter=',')
            kD    = np.genfromtxt(path_kD, delimiter=',')
            kW    = np.genfromtxt(path_kW, delimiter=',')

            print('done')
            simtD = funPS_D_t(freqt)
            simtW = funPS_W_t(freqt)
            np.save(f"/rds/general/user/csheeran/home/sim_dat/simtD_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.npy", simtD)
            np.save(f"/rds/general/user/csheeran/home/sim_dat/simtW_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.npy", simtW)
            print(simtD)
            print(simtW)
            print(tD)
            print('nearby')

            time = plt.figure(1)
            plt.scatter(freqt, tD, label='Sim Drive')
            plt.scatter(freqt, tW, label='Sim Wild')
            plt.scatter(freqt, simtD, label='MEQ Drive')
            plt.scatter(freqt, simtW, label='MEQ Wild')
            plt.xscale('log')
            plt.yscale('log')
            plt.legend()
            plt.title('Logscale Plot with Scatter Plots')
            plt.xlabel('Temporal Frequency')
            plt.ylabel('Power')
            time.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_time_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
            time.clf()            
            
            #spacetime = plt.figure(2)
            #plt.scatter(freqk, kD, label='Sim Drive')
            #plt.scatter(freqk, kW, label='Sim Wild')
            #plt.scatter(freqk, simkD, label='MEQ Drive')
            #plt.scatter(freqk, simkW, label='MEQ Wild')
            #plt.xscale('log')
            #plt.yscale('log')
            #plt.legend()
            #plt.title('Logscale Plot with Scatter Plots')
            #plt.xlabel('Spatial Frequency')
            #plt.ylabel('Power')
            #spacetime.savefig(f"/rds/general/user/csheeran/home/PS_plots/PS_space_h={round_h}_s={round_s}_ep={round_ep}_Cc={round_Cc}.pdf", bbox_inches='tight')
