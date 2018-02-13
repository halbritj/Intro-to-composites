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
