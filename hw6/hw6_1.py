import numpy as np
import matplotlib.pyplot as plt



n_0 = 7.93 * 10**-14 #Pa*S
U = 9.08 * 10**4 #J/mol
R = 8.3143 #J/mol*K
K = 14.1

n = np.vectorize(lambda T, alpha: n_0*np.exp(U/(R*T) + K*alpha))

T = np.array([350,375,400], float)
alpha = np.linspace(0,.5,100)

n_T = n(T[:, None], alpha)

style = ['k--', 'k:', 'k']


fig = plt.figure()
fig.suptitle(r'Viscocity vs degree of cure')

ax = fig.add_subplot(111)

ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$\eta(T,\alpha)$  ($Pa.S$)') 


for i in range(3):
    label = 'T=%d K' %T[i]
    ax.plot(alpha, n_T[i], style[i], label=label)


ax.legend()
plt.show()

"""
The graph reveals the inverse relationship between temperature and viscocity.
Higher temperature cures result in a matrix with a much lower viscocity of the matrix.
Here a difference of 50K between 350K and 400K results in a matrix that is
50x less viscous. The lower viscocity is caused by the temperatue of matrix
weaking the intermolecular bonds (vaan der waals/hydrogen bonding).
"""
