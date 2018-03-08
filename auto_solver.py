import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=2)


def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

Em = .4*10**6 #psi
Ef = 30*10**6 #psi

vf12 = .3
vm = .5

h = .15*10**-3

E1_ = lambda Vf: Ef*Vf + Em*(1-Vf)

def E2_(Vf):
    a = np.sqrt(Vf)
    return ( (1-a)/Em + a/(Ef*a + Em*(1-a)) )**-1

v12_ = lambda Vf: vf12*Vf + vm*(1-Vf)

theta = np.array([0,45,-45,-45,45,0])

x = []
y = []

for Vf in np.linspace(.6, .8, 100):
    E1 = E1_(Vf)
    E2 = E2_(Vf)    
    v12 = v12_(Vf)
    G12 = E1/(2*(1-v12))
    
    S = np.array([
        [1/E1, -v12/E1, 0],
        [-v12/E1, 1/E2, 0],
        [0,  0       , 1/G12]], np.float64)



    N = theta.size

    H = N*h

    Z = np.arange(N+1)*h - .5*H
    
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


    x.append(Vf)
    y.append(-abd[0,1]/abd[0,0])

#plt.plot(x,y)
#plt.show()
    
