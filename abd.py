import numpy as np
import sympy as sp
#np.set_printoptions(precision=4)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

theta = np.array([90,45,0,-45])
    
mod = [2, 's']
for m in mod:
    if isinstance(m, int):
        for i in range(m-1):
            theta = np.hstack((theta, theta))
    if isinstance(m, str):
        if m == 's':
            theta = np.hstack((theta, theta[::-1]))

#theta = np.array([0,60,-60,-60,60,0])




def getABD(theta, h):

    N = theta.size

    H = N*h

    Z = np.arange(N+1)*h - .5*H

    E1 = 155 * 10**9
    E2 = 12.1 * 10**9
    v12 = .2
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

    return ABD, np.linalg.inv(ABD)

h = .15*10**-3
ABD, abd = getABD(theta, h)

E1, E2, v12, v21, G12 = sp.var('E1, E2, v12, v21, G12')


def equ(rhs, lhs):
    return rhs - lhs

equs = [
    E1/v12 - E2/v21,
    G12 - E1/(2*(1-v12))
    ]

subs = {
    E1: 1,
    E2: 2,
    v12: 3,
    }

def evalf(f, subs):
    for var in subs:
        value = subs[var]
        f = f.subs(var, value)
    return f

equs = [evalf(equ, subs) for equ in equs]




