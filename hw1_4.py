import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


def L(plane):
    if plane == '12': return sp.diag(1,1,-1)
    if plane == '23': return sp.diag(-1,1,1)
    if plane == '31': return sp.diag(1,-1,1)
    
s1, s2, s3, s4, s5, s6 = sp.var('s1 s2 s3 s4 s5 s6')

Stress = sp.Matrix([
    [s1, s6, s5],
    [s6, s2, s4],
    [s5, s4, s3]])

C = sp.Matrix([
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0],
    [0,0,0,0,0,0]])
    
C_vars = []
for i in range(6):
    for j in range(i+1):
        c = sp.var('C%d%d' %(j+1, i+1))
        C_vars.append(c)
        C[i,j] = c
        C[j,i] = c


S = sp.Matrix([
    [0],
    [0],
    [0],
    [0],
    [0],
    [0]])
S_vars = []
for i in range(6):
    s = sp.var('S%d' %(i+1))
    S_vars.append(s)
    S[i,0] = s


E = sp.Matrix([
    [0],
    [0],
    [0],
    [0],
    [0],
    [0]])
E_vars = []
for i in range(6):
    e = sp.var('E%d' %(i+1))
    E_vars.append(e)
    E[i,0] = e



    
Sp = sp.Matrix([
    [0],
    [0],
    [0],
    [0],
    [0],
    [0]])
Sp_vars = []
for i in range(6):
    s = sp.var('Sp%d' %(i+1))
    Sp_vars.append(s)
    Sp[i,0] = s


Ep = sp.Matrix([
    [0],
    [0],
    [0],
    [0],
    [0],
    [0]])
Ep_vars = []
for i in range(6):
    ep = sp.var('Ep%d' %(i+1))
    Ep_vars.append(ep)
    Ep[i,0] = ep


e = sp.diag(1,1,1,-1,-1,1)


for i in range(6):
    Ep = Ep.subs(Ep_vars[i], e[i,i]*E_vars[i])

out = S
S = C*E
Sp = C*Ep
for i in range(6):
    out[i] = S[i] - e[i,i]*Sp[i]

