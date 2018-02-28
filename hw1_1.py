import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

C = np.array([
    [ 158.   ,  5.64 ,  5.64 ,  0.   ,  0.   ,  0.  ],
    [   5.64 , 15.51 ,  7.21 ,  0.   ,  0.   ,  0.  ],
    [   5.64 ,  7.21 , 15.51 ,  0.   ,  0.   ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  3.2  ,  0.   ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  0.   ,  4.4  ,  0.  ],
    [   0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  4.4 ]], np.float64) * 10**9

def Transform(theta):
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

theta = np.linspace(-90, 90, 100)

T = Transform(theta)
T_ = np.rollaxis(T, 2)

#S_bar = np.einsum('...jk,kl,...lm->...jm', T.T, s, T_) #[S_bar] = [T.T][S][T]
S_bar = np.matmul(T.T, np.matmul(s, T_))


Q_bar = np.linalg.inv(S_bar)
    
eps_xyz = np.array([[0,0,1]]).T

sigma_xyz = np.dot(Q_bar, eps_xyz).reshape(100, -1)

Q_bar /= 10**9

fig, ax = plt.subplots()

plt.plot(theta, Q_bar[:,0,2], 'k-', label=r'$\bar Q_{16}$')
plt.plot(theta, Q_bar[:,1,2], 'k--', label=r'$\bar Q_{26}$')
plt.plot(theta, Q_bar[:,2,2], 'k:', label=r'$\bar Q_{66}$')

plt.title(r'Stiffness factors relevant to $\gamma_{xy}$')
plt.xlabel(r'$\theta^\circ$', fontsize=15)
plt.ylabel(r'$GPa$', fontsize=15)

legend = ax.legend(loc='upper right', shadow=True)
plt.xticks(np.linspace(-90, 90, 13))

plt.show()
