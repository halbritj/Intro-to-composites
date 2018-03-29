import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)

def Transform(theta):
    m = np.cos( np.deg2rad(theta) )
    n = np.sin( np.deg2rad(theta) )
    return np.array([
        [m**2, n**2, 2*m*n],
        [n**2, m**2, -2*m*n],
        [-m*n, m*n,  m**2 - n**2]], np.float64)


theta = np.array([+30,-30,0,0,-30,+30])

N = theta.size
h = .15*10**-3
H = N*h

Z = np.arange(N+1)*h - .5*H

E1 = 155 * 10**9
E2 = 12.1 * 10**9
v12 = .248
G12 = 4.4 * 10**9

S = np.array([
    [1/E1, -v12/E1, 0],
    [-v12/E1, 1/E2, 0],
    [0,  0       , 1/G12]], np.float64)

T = Transform(theta)
T_ = np.rollaxis(T, 2)

Sbar = np.einsum('...jk,kl,...lm->...jm', T.T, S, T_)
Qbar = np.linalg.inv(Sbar)

A = np.sum( np.diff(Z)[:, None, None] * Qbar, axis=0)
B = (1/2)*np.sum( np.diff(Z**2)[:, None, None] * Qbar, axis=0)
D = (1/3)*np.sum( np.diff(Z**3)[:, None, None] * Qbar, axis=0)

ABD = np.block([[A,B],[B,D]])

ABD[np.abs(ABD) < 10**-8] = 0

abd = np.linalg.inv(ABD)

Xt = 1500 * 10**6
Xc = 1250 * 10**6

Yt =   50 * 10**6
Yc =  200 * 10**6
S  =  100 * 10**6


alpha_1 = -0.018 * 10**-6
alpha_2 = 24.3   * 10**-6

beta_1 =  146 * 10**-6
beta_2 = 4770 * 10**-6


alpha = np.array([alpha_1, alpha_2, 0])
beta  = np.array([beta_1, beta_2, 0])

alpha_bar = np.matmul(T.T, alpha[:, None]).reshape(N, 1, 3)
beta_bar = np.matmul(T.T, beta[:, None]).reshape(N, 1, 3)

NM = np.array([[0,0,0,1,0,0]]).T

N_t = np.sum( np.diff(Z)[:, None] *  np.sum( Qbar * alpha_bar, axis=2 ), axis=0)
M_t = (1/2)*np.sum( np.diff(Z**2)[:, None] *  np.sum( Qbar * alpha_bar, axis=2 ), axis=0)
NM_t = np.concatenate([N_t, M_t]).reshape(6, 1)

N_m = np.sum( np.diff(Z)[:, None] *  np.sum( Qbar * beta_bar, axis=2 ), axis=0)
M_m = (1/2)*np.sum( np.diff(Z**2)[:, None] *  np.sum( Qbar * beta_bar, axis=2 ), axis=0)

NM_m = np.concatenate([N_m, M_m]).reshape(6, 1)

dT = -150

NM_r = 25*NM + dT*NM_t

Z_ = np.hstack((Z[:N//2], Z[-N//2:]))

ref_surf = np.dot(abd, NM_r)
surf_strain = ref_surf[:3]
surf_curve = ref_surf[3:]

Strains = surf_strain + np.dot(surf_curve, Z_[None,:])
Sigma_xyz = np.matmul(Qbar, Strains.T[:,:,None])
Sigma_123 = np.matmul(T_, Sigma_xyz)


Z_stack = np.array([np.linspace(i*h, (i+1)*h) - .5*H for i in range(N)])
strains = surf_strain[None, :] + np.einsum('ij,kl->ikj', Z_stack, surf_curve)

sigma_xyz = np.einsum('ijk,ikl->ijl', Qbar, strains)
sigma_123 = np.einsum('ijk,ikl->ijl', T_, sigma_xyz)

a = np.matmul( T_, np.matmul(Qbar, surf_curve) )


Z_stack *= 10**3
Z_ *= 10**3
Z *= 10**3
H *= 10**3



f, ax = plt.subplots(1, 3, sharey=True)

ax[0].set_ylim([-H/2, H/2])

ax[0].set_title(r'$\sigma_1$')
ax[0].axvline(-Xc/10**6, dashes=(1,2), color = (1,0,0))
ax[0].axvline(Xt/10**6, dashes=(1,2), color = (1,0,0))
ax[0].text(-Xc/10**6, H/2, r'$X_c$')
ax[0].text(Xt/10**6, H/2, r'$X_t$')

ax[1].set_title(r'$\sigma_2$')
ax[1].axvline(-Yc/10**6, dashes=(1,2), color = (1,0,0))
ax[1].axvline(Yt/10**6, dashes=(1,2), color = (1,0,0))
ax[1].text(-Yc/10**6, H/2, r'$Y_c$')
ax[1].text(Yt/10**6, H/2, r'$Y_t$')

ax[2].set_title(r'$\tau_{12}$')
ax[2].axvline(-S/10**6, dashes=(1,2), color = (1,0,0))
ax[2].axvline(S/10**6, dashes=(1,2), color = (1,0,0))
ax[2].text(-S/10**6, H/2, r'$-S$')
ax[2].text(S/10**6, H/2, r'$S$')


for i in range(3):
    ax[i].axvline(0, dashes=(1,2))
    for z in Z:
        ax[i].axhline(z, dashes=(1,2))

for i in range(N):
    for j in range(3):
        ax[j].plot(sigma_123[i,j,:] / 10**6, Z_stack[i])

        ax[j].scatter( Sigma_123[i,j,0] / 10**6, Z_[i])

s1 = sigma_123[:, 0]
s2 = sigma_123[:, 1]
s3 = sigma_123[:, 2]

axial = np.abs(s1) / ( (s1 >= 0)*Xt + (s1 < 0)*Xc )
trans = np.abs(s2) / ( (s2 >= 0)*Yt + (s2 < 0)*Yc )
shear = np.abs(s3) / S
        
plt.show()

'''
>>> Sigma_123 / 10**6
array([[[-210.46239517],
        [  -1.9220824 ],
        [  24.54318073]],

       [[-291.11679885],
        [  -1.28690287],
        [  -8.46186138]],

       [[-145.04600046],
        [ -14.22428134],
        [   2.14382016]],

       [[ 230.75092989],
        [ -45.10693709],
        [  -2.14382016]],

       [[ 159.62084854],
        [ -45.04191657],
        [  30.13757409]],

       [[  78.96644486],
        [ -44.40673704],
        [ -46.21889343]]])
'''
