import numpy as np

np.set_printoptions(precision=3)

sqrt = lambda e: np.sqrt(e)

Ef1 = 80*10**9
Ef2 = 80*10**9
Gf12 = 33.33*10**9
vf12 = .2
Vf = .6
sq_Vf = sqrt(Vf)

Em = 3.35*10**9
Gm = 1.24*10**9
vm = .35
Vm = 1 - Vf


E1 = (Ef1*Vf + Em*Vm)
E2 = ( (1-sq_Vf)/Em + sq_Vf/(Ef2*sq_Vf + Em*(1-sq_Vf)) )**-1
G12 = Gm * ( (Gm*(1-Vf) + Gf12*(1+Vf)) / (Gm*(1+Vf) + Gf12*(1-Vf)) )
v12 = vf12*Vf + vm*Vm


S = np.array([
    [1/E1       , -v12/E1   , 0],
    [-v12/E1    , 1/E2      , 0],
    [0          , 0         , 1/G12]], np.float64)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)

T = Transform(0)
Sbar_0 = np.dot(T.T, np.dot(S, T))

T = Transform(45)
Sbar_45 = np.dot(T.T, np.dot(S, T))
