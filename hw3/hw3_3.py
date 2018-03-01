import numpy as np

#np.set_printoptions(precision=2)

m = lambda a: np.cos(np.deg2rad(a))
n = lambda a: np.sin(np.deg2rad(a))

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

theta = np.array([0,90,90,0])

N = theta.size
h = 150*10**-6
H = N*h

Z = np.arange(N+1)*h - .5*H

alpha1 = -.018*10**-6
alpha2 = 24.3*10**-6
alpha3 = alpha2

alpha = np.array([[alpha1, alpha2, 0]]).T

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

#alpha_bar = np.matmul(T.T, alpha)
alpha_bar = np.array([2.1, 2.1, 0])*10**-6


Sbar = np.einsum('...jk,kl,...lm->...jm', T.T, S, T_)
Qbar = np.linalg.inv(Sbar)

#N_t = np.sum( np.diff(Z)[:, None, None] * np.matmul(Qbar, alpha_bar), axis=0 )
#M_t = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * np.matmul(Qbar, alpha_bar), axis=0 )

N_t = np.sum( np.diff(Z)[:, None] * np.dot(Qbar, alpha_bar), axis=0 )
M_t = (1/2)*np.sum( np.diff(Z**2)[:, None] * np.dot(Qbar, alpha_bar), axis=0 )



Sbar = np.einsum('...jk,kl,...lm->...jm', T.T, S, T_)
Qbar = np.linalg.inv(Sbar)

A = np.sum( np.diff(Z)[:, None, None] * Qbar, axis=0)
B = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * Qbar, axis=0)
D = (1/3)*np.sum( np.diff(Z**3)[:, None, None] * Qbar, axis=0)

ABD = np.vstack( (np.hstack((A,B)), np.hstack((B,D))) )
ABD[np.abs(ABD) < 10**-8] = 0

abd = np.linalg.inv(ABD)

strains = np.dot(abd[:2, :2], N_t[:2])

