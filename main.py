import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


C = np.array([
    [ 158.   ,  5.64 ,  5.64 ,  0.   ,  0.   ,  0.  ],
    [   5.64 , 15.51 ,  7.21 ,  0.   ,  0.   ,  0.  ],
    [   5.64 ,  7.21 , 15.51 ,  0.   ,  0.   ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  3.2  ,  0.   ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  0.   ,  4.4  ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  4.4 ]], np.float64) * 10**9


def T(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

S = np.linalg.inv(C)

s = np.zeros((3,3), np.float64)
s[:2,:2] = S[:2, :2]
s[-1, -1] = S[-1, -1]

#Q = np.linalg.inv(S_plane)


theta = np.linspace(-90, 90, 100)


T_ = np.rollaxis(T(theta), 2)

s_bar = np.ndarray((100,3,3))
for i, t in enumerate(T_):
    s_bar[i] = np.dot(t.T, np.dot(s, t))

Q_bar = np.linalg.inv(s_bar)
    
eps_xyz = np.array([[0,0,1]]).T

sigma_xyz = np.dot(Q_bar, eps_xyz).reshape(100, -1)

plt.plot(theta, sigma_xyz[:, 0]) #blue
plt.plot(theta, sigma_xyz[:, 1]) #green
plt.plot(theta, sigma_xyz[:, 2]) #red
plt.show()





'''
E1, E2, E3 = S.diagonal()[:3]**-1
V12 = -S[0,1] * E1
V21 = (E2/E1) * V12

V23 = -S[1,2] * E2
V32 = (E3/E2) * V23

V13 = -S[0,2] * E1
V31 = (E3/E1) * V13

G12= S[5,5]**-1

delta = (1 - V12*V21 - V23*V32 - V13*V31 - 2*V21*V32*V13)/(E1*E2*E3)


E1,E2,E3,V12,V13,V23 = sp.var('E1 E2 E3 V12 V13 V23')

E = sp.diag(1/E1, 1/E2, 1/E3)

S = E * sp.Matrix([
                [1, -V12, -V13],
                [-V12, 1, -V23],
                [-V13, -V23, 1]])
'''

