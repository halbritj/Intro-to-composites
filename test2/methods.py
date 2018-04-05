import numpy as np

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

def getLaminate(theta, h, E1, E2, v12, G12):
    N = theta.size
    H = N*h

    Z = np.arange(N+1)*h - .5*H

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

    return ABD, abd, Sbar, Qbar, T, T_









