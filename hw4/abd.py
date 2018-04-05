import numpy as np
import sympy as sp

np.set_printoptions(precision=4)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)



#theta = np.array([0,90,45,-45,-45,45,90,0])
#theta = np.array([0,45,90,-45,0,45,90,-45,-45,90,45,0,-45,90,45,0])
theta = np.array([15,-15,15,-15,15,-15,15,-15,15,-15,15,-15,15,-15,15,-15])

N = theta.size
h = .15*10**-3
H = N*h

Z = np.arange(N+1)*h - .5*H

E1 = 155 * 10**9
E2 = 12.1 * 10**9
v12 = .248
G12 = 4.4 * 10**9

S = np.array([
    [1/E1, -v12/E1, 0],
    [-v12/E1, 1/E2, 0],
    [0,  0       , 1/G12]], np.float64)

T = Transform(theta)
T_ = np.rollaxis(T, 2)

Sbar = np.einsum('...jk,kl,...lm->...jm', T.T, S, T_)
Qbar = np.linalg.inv(Sbar)

A = np.sum( np.diff(Z)[:, None, None] * Qbar, axis=0)
B = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * Qbar, axis=0)
D = (1/3)*np.sum( np.diff(Z**3)[:, None, None] * Qbar, axis=0)

ABD = np.vstack( (np.hstack((A,B)), np.hstack((B,D))) )
ABD[np.abs(ABD) < 10**-8] = 0

abd = np.linalg.inv(ABD)
