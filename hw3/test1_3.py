import numpy as np
import sympy as sp

np.set_printoptions(precision=2)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)


theta = np.array([0,90,53,-53])
theta = np.hstack((theta, theta))
theta = np.hstack((theta, theta[::-1]))

E1 = 155 *10**9 #GPa -> Pa
E2 = 12.1 *10**9 #GPa -> Pa
v12 = .248
G12 = 4.40 *10**9 #GPa -> Pa

h = .15 *10**-3 #mm -> m

alpha1 = -.018 *10**-6 #1/C
alpha2 = 24.3 *10**-6 #1/C

beta1 = 146 *10**-6 #1/%m
beta2 = 4770 *10**-6 #1/%m

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

alpha = np.array([alpha1, alpha2, 0])
alpha_bar = np.matmul(T.T, alpha[:, None])

beta = np.array([beta1, beta2, 0])
beta_bar = np.matmul(T.T, beta[:, None])

N_t = np.sum( np.diff(Z)[:, None, None] * np.matmul(Qbar, alpha_bar), axis=0 )
M_t = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * np.matmul(Qbar, alpha_bar), axis=0 )

N_m = np.sum( np.diff(Z)[:, None, None] * np.matmul(Qbar, beta_bar), axis=0 )
M_m = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * np.matmul(Qbar, beta_bar), axis=0 )


TR = np.vstack((N_t, M_t))
MR = np.vstack((N_m, M_m))

NR = np.array([[.2*10**6, -.2*10**6, 0, 0, 0, 0]]).T


dM = .5


i = 1
def sigma(dT):
    N = NR + TR*dT + MR*dM
    a = np.dot(abd, N)[:3]
    #print(a)
    temp = a - alpha_bar[i]*dT - beta_bar[i]*dM
    #print(temp)
    result = np.dot(T_[i], np.dot(Qbar[i], temp))
    print(result)
    


a = sigma(-125)



