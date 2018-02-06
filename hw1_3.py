import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

np.set_printoptions(precision=3)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)


#solve problem assuming plane stress

E1 = 50 * 10**9
E2 = 15.2 * 10**9
v12 = 0.254
G12 = 4.70 *10**9

#V12/E1 = V21/E2


S = np.array([
    [1/E1       , -v12/E1   , 0],
    [-v12/E1    , 1/E2      , 0],
    [0          , 0         , 1/G12]], np.float64)

theta = np.linspace(-90, 90, 100)

t = Transform(theta)
T_ = np.rollaxis(t, 2)


T = t.T


S_bar = np.dot(T_, S)


S11_bar = S_bar[:,0,0]
S22_bar = S_bar[:,1,1]
S16_bar = S_bar[:,0,2]
S26_bar = S_bar[:,2,1]

nu_xy_x = S16_bar / S11_bar
nu_xy_y = S26_bar / S22_bar

plt.plot(theta, S11_bar)
#plt.plot(theta, S22_bar)

#plt.plot(theta, S16_bar)
#plt.plot(theta, S26_bar)

#plt.plot(theta, nu_xy_x)
#plt.plot(theta, nu_xy_y)
#plt.show()

'''
E1, E2, V12, G12 = sp.var('E1 E2 V12 G12')
m, n = sp.var('m n')


T = sp.Matrix([
    [],
    [],
    []])

S = sp.Matrix([
    [1/E1, -V12/E1, 0],
    [-V12/E1, 1/E2, 0],
    [0, 0, 1/G12]])
'''



