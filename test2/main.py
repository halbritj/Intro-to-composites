import numpy as np
from methods import *
from sympy import *


size = 3
upper = np.array([[0]*(i) + [1]*(size-i) for i in range(size)], bool)

theta = np.array([15,-15,15,-15,15,-15,15,-15,15,-15,15,-15,15,-15,15,-15])

h = .15*10**-3
E1 = 155 * 10**9
E2 = 12.1 * 10**9
v12 = .248
G12 = 4.4 * 10**9

ABD, abd, Sbar, Qbar, T, T_ = getLaminate(theta, h, E1, E2, v12, G12)

A = ABD[:3, :3]
B = ABD[-3:, :3]
D = ABD[-3:, -3:]

A11, A12, A16, A22, A26, A66 = A[upper]
B11, B12, B16, B22, B26, B66 = B[upper]
D11, D12, D16, D22, D26, D66 = D[upper]

x,y = var('x y')
a,b = var('a b')
m,n = var('m n')


q0 = (4500*x*(a-x)**3)/(a*b) #N/m/m

equ = q0 * sin(m*pi*x/a) * sin(n*pi*y/b)
Q = integrate( integrate(equ, (x,0,a)), (y,0,b)) * 4 * (a*b)**-1


a = 1 #m
b = 1 #m

subs = {
    'a': a,
    'b': b}

Q = Q.subs(subs)


def Qmn(m, n): #uniform loading
    result = np.zeros(m.shape)

    size, size = m.shape

    for i in range(size):
        for j in range(size):
            result[i,j] = Q.evalf(subs={'m': m[i,j], 'n': n[i,j]})

    #return float(Q.evalf(subs=subs))
    return result
    
    #return (16*q0)/(np.pi**2 * m * n)

def mn(m_range, n_range):
    m, n = np.meshgrid( np.arange(m_range)+1, np.arange(n_range)+1 )
    return m, n

def a_mn(m, n, Nx=0, Ny=0): #Navier SS2 case
    A = (m*np.pi)/a
    B = (n*np.pi)/b
    
    c11 =  A11*(A**2) + A66*(B**2)
    c12 = (A12 + A66)*A*B
    c13 = -(3*B16*(A**2) + B26*(B**2))*B
    c22 = A66*(A**2) + A22*(B**2)
    c23 = -(B16*(A**2) + 3*B26*(B**2))*A
    c33 = D11*(A**4) + 2*(D12 + 2*D66)*(A**2)*(B**2) + D22*(B**4)

    a0 = c11*c22 - c12*c12
    a1 = c12*c23 - c13*c22
    a2 = c13*c12 - c11*c23

    amn = c33 + (c13*a1 + c23*a2)/a0

    return amn

def Wmn(m, n):
    return Qmn(m, n) / a_mn(m, n)


def w_o(x, y, precision=3):
    W = 0
    W_new = 1

    size = 1
    while size < 15:
        p = np.floor( -np.log10(np.abs(W_new-W)) )
        if p == np.inf: p = precision
        
        print('size: %d \tprecision: %d\tw: %.10f' %(size, p, W_new))

        W = W_new
        
        m, n = mn(size, size)

        #print(Wmn(m,n))

        W_new = np.sum( Wmn(m, n) * np.sin(m*np.pi*x/a) * np.sin(n*np.pi*y/b) )
    
        size += 1

    return W_new

np.seterr(divide='ignore')



w = []

for x in np.linspace(.4, .5, 20):
    r = w_o(x, .4, precision=5)
    w.append(r)

import matplotlib.pyplot as plt

plt.plot(np.linspace(.4, .5, 20), w)
plt.show()


