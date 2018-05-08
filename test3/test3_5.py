import numpy as np
import matplotlib.pyplot as plt

class resin:
    def __init__(self, name, params):
        self.name = name

        self.A_1, self.E_A1, self.m_1, self.n_1, self.A_2, self.E_A2, \
        self.m_2, self.n_2, self.D, self.alpha_C0, self.alpha_CT = params



A = resin('MTM45-1', [
    25.3*10**4,
    60628,
    0.55,
    21.11,
    4.84*10**4,
    61752,
    0.80,
    1.18,
    44.3,
    -1.4,
    5.33*10**-3])
           
B = resin('5320', [
    8.23*10**7,
    82375,
    0.75,
    12.46,
    1.04*10**5,
    62355,
    0.90,
    2.07,
    40.4,
    -1.12,
    4.53*10**-3])

T = 413 #K
R = 8.31432 #J/mol.K
def cure_rate(alpha, r):
    K1 = r.A_1 * np.exp(-r.E_A1 / (R*T))
    K2 = r.A_2 * np.exp(-r.E_A2 / (R*T))

    a = K1*(alpha**r.m_1)*((1-alpha)**r.n_1)
    b = K2*(alpha**r.m_2)*((1-alpha)**r.n_2)
    c = 1 + np.exp(r.D*(alpha-(r.alpha_C0 + r.alpha_CT*T)))


    return a + b/c


alpha = np.linspace(0,1,1000)

A_cure = cure_rate(alpha, A)
B_cure = cure_rate(alpha, B)



fig = plt.figure()
fig.suptitle('Cure Rate vs Degree of Cure  (T=413K)')
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$d\alpha$/$dt$ $[s^{-1}]$') 

ax.plot(alpha, A_cure, 'k--', label='MTM45-1')
ax.plot(alpha, B_cure, 'k:', label='5320')

ax.axvline(.3, ls='--', alpha=.3)
ax.axvline(.6, ls='--', alpha=.3)

ax.legend(title='Epoxy Resin\n   System')
plt.show()




