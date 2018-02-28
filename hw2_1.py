import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

Erandom = (1+.5)*1.25*10**6 #psi

Wf = .2
rho_f = 1.8 #g/cm**3
Ef = 30*10**6 #psi
lf = sp.var('l_f')
df = .0006 #in

Wm = 1 - Wf
rho_m = 1.14 #g/cm**3
Em = .4*10**6 #psi

rho_c = 1/( (Wf/rho_f) + (Wm/rho_m) )
Vf = Wf * (rho_c/rho_f)

eta_L = ((Ef/Em)-1) / ( (Ef/Em) + 2*(lf/df))
eta_T = ((Ef/Em)-1) / ( (Ef/Em) + 2)

E1 = ( (1 + 2*(lf/df)*eta_L*Vf) / (1 - eta_L*Vf) )*Em
E2 = ( (1+2*eta_T*Vf) / (1-eta_T*Vf) )*Em

f = .375*E1 + .625*E2 - Erandom
length, = sp.solve(f, lf)
E = sp.lambdify(lf, f+Erandom)
lf = np.linspace(0, 1, 1000)
fig, ax = plt.subplots()

pair = (length, E(length)/10**6)
ax.plot(*pair, 'o')
ax.annotate(r'$Min$ $\l_f={%.3f}$ $in$' %pair[0], xy=(length, pair[1]-.1), fontsize=13)
plt.plot(lf, E(lf)/10**6, label='$E_{random}(l_f)$')

plt.axhline(Erandom/10**6, linestyle='--')
ax.annotate(r'$E_{random}=%.3f*10^6$ $psi$' %(Erandom/10**6), xy=(.6, Erandom/10**6 - .07), fontsize=13)

plt.title('$E_{random}$ vs fiber length')
plt.xlabel(r'$inches$', fontsize=15)
plt.ylabel(r'$10^6$ $psi$', fontsize=15)

legend = ax.legend(loc='lower right', shadow=True)

plt.show()
