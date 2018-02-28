import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
#variables relevant to equation
E1, E2, m, n, v12, G12, R = sp.var('E1 E2 m n v12 G12 R')

#shear modulus equation
G12 = (R*E1)/(2*(1+v12))

#Transformed tensile modulus
Ex = E1/(m**4 + (n**2)*(m**2)*(E1/G12 - 2*v12) + (E1/E2)*(n**4))

#normalize modulus, equation to be plotted as function of theta
f = Ex/E1

#assume volume fractions
Vf = .625
Vm = 1 - Vf

#assume fiber poisson's ratios
vf12 = .215
vm = .325
v12 = vf12*Vf + vm*Vm
E1E2 = 13.8

theta = np.linspace(0,90,90)
m = np.cos(np.deg2rad(theta))
n = np.sin(np.deg2rad(theta))

fig, ax = plt.subplots()

R = .75
ExE1_1 = 1/(E1E2*n**4 + m**4 + m**2*n**2*(-2*v12 + (2*v12 + 2)/R))
plt.plot(theta, ExE1_1, 'k-', label=r'$R=.75$')

R = 1.25
ExE1_2 = 1/(E1E2*n**4 + m**4 + m**2*n**2*(-2*v12 + (2*v12 + 2)/R))
plt.plot(theta, ExE1_2, 'k--', label=r'$R=1.25$')

plt.axhline(1/E1E2, linestyle='-.')
ax.annotate(r'${(E_1/E_2)}^{-1}=%.3f$' %(1/13.8), xy=(10, 1/12.5), fontsize=13)

plt.title('Normalized Modulus vs Fiber Orientation Angle\n(Multiple $R$ values)')
plt.xlabel(r'$\theta^\circ$', fontsize=15)
plt.ylabel(r'$E_x/E_1$', fontsize=15)
legend = ax.legend(loc='upper right', shadow=True)
plt.xticks(np.linspace(0, 90, 10))
plt.show()
