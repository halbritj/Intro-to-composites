import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

m, n, t = sp.var('m n t')

T_inv = sp.Matrix([
    [m**2, n**2, -2*m*n],
    [n**2, m**2, 2*m*n],
    [m*n, -m*n, m**2-n**2]])

'''
subs={m: sp.cos(t), n: sp.sin(t)}
for var in subs:
    T_inv = T_inv.subs(var, subs[var])
'''
Q11, Q12, Q22, Q66 = sp.var('Q11, Q12, Q22, Q66')

Q = sp.Matrix([
    [Q11, Q12, 0],
    [Q12, Q22, 0],
    [0  , 0  , Q66]])

Q_bar = (T_inv * Q * T_inv.T)



ply = np.array([90,45,0,-45])
    
mod = [2, 's']
for m in mod:
    if isinstance(m, int):
        for i in range(m-1):
            ply = np.hstack((ply, ply))
    if isinstance(m, str):
        if m == 's':
            ply = np.hstack((ply, ply[::-1]))

F = Q_bar[1,2]
m, n, t = sp.var('m n t')

subs = {m: sp.cos(t), n: sp.sin(t)}
for key in subs:
    F = F.subs(key, subs[key])

theta = np.array([0,0,22.5,45,-45,-22.5,90,90])

h = .15*10**-3
N = 2*theta.size
H = N*h

Z = np.arange(N+1)*h - .5*H

z3 = np.diff(Z**3)
#plt.plot(z3)
#plt.show()

k,h,N = sp.var('k,h,N')

f = ( (k+1)*h - .5*N*h )**3 - (k*h - .5*N*h)**3



