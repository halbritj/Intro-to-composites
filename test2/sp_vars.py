from sympy import *
from methods import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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



N = theta.size
H = N*h

Z = np.arange(N+1)*h - .5*H
Z_ = np.hstack((Z[:N//2], Z[-N//2:]))




x,y = var('x y')
m,n = var('m n')

a,b = 1,1

q0 = (4500*x*(a-x)**3)/(a*b) #N/m/m
equ = q0 * sin(m*pi*x/a) * sin(n*pi*y/b)
Qmn = integrate( integrate(equ, (x,0,a)), (y,0,b)) * 4 / (a*b)

Qmn = lambdify((x,y,m,n), Qmn)

C_m = lambdify((x, m), cos(m*pi*x/a))
S_m = lambdify((x, m), sin(m*pi*x/a))

C_n = lambdify((y, n), cos(n*pi*y/b))
S_n = lambdify((y, n), sin(n*pi*y/b))

u_trig = cos(m*pi*x/a)*sin(n*pi*y/b)
v_trig = sin(m*pi*x/a)*cos(n*pi*y/b)
w_trig = sin(m*pi*x/a)*sin(n*pi*y/b)

symbols = (x,y,m,n)


Strain = lambdify(symbols, Matrix([
    [diff(u_trig, x), diff(u_trig, y)],
    [diff(v_trig, x), diff(v_trig, y)]]))

Curve = lambdify(symbols, Matrix([
    [diff(w_trig, x, x), diff(w_trig, x, y)],
    [diff(w_trig, y, x), diff(w_trig, y, y)]]))

u_trig = lambdify(symbols, u_trig)
v_trig = lambdify(symbols, v_trig)
w_trig = lambdify(symbols, w_trig)




def evals(x, y, m, n):
    Qmn_ = Qmn(x, y, m, n)
    
    A = m*np.pi/a
    B = n*np.pi/b

    c11 =  A11*(A**2) + A66*(B**2)
    c12 = (A12 + A66)*A*B
    c13 = -(3*B16*(A**2) + B26*(B**2))*B
    c22 = A66*(A**2) + A22*(B**2)
    c23 = -(B16*(A**2) + 3*B26*(B**2))*A
    c33 = D11*(A**4) + 2*(D12 + 2*D66)*(A**2)*(B**2) + D22*(B**4)

    a0 = c11*c22 - c12*c12
    a1 = c12*c23 - c13*c22
    a2 = c13*c12 - c11*c23

    a_mn = c33 + (c13*a1 + c23*a2)/a0

    Wmn = Qmn_ / a_mn
    Umn = (a1/a0)*Wmn
    Vmn = (a2/a0)*Wmn

    return Wmn, Umn, Vmn

def disp(x, y, size):
    #u = Umn * C_m * S_n
    #v = Vmn * S_m * C_n
    #w = Wmn * S_m * S_n

    u,v,w = 0,0,0    
    for m in range(1, size+1):
        for n in range(1, size+1):
            args = (x,y,m,n)
            Wmn, Umn, Vmn = evals(*args)
            u += Umn * u_trig(*args)
            v += Vmn * v_trig(*args)
            w += Wmn * w_trig(*args)
            
    return u.evalf(), v.evalf(), w.evalf()


def strains(x, y, size):
    E_x, E_y, E_xy = 0, 0, 0
    K_x, K_y, K_xy = 0, 0, 0

    S_total = np.zeros((2,2))
    C_total = np.zeros((2,2))
    

    for m in range(1, size+1):
        for n in range(1, size+1):
            args = (x,y,m,n)
            Wmn, Umn, Vmn = evals(*args)

            S_total += Strain(*args) * [[Umn],[Vmn]]
            C_total += Curve(*args) * Wmn


    (u_x, u_y), (v_x, v_y) = S_total
    (w_xx, w_xy), (_, w_yy) = C_total
    
    return u_x, v_y, v_x+u_y, -w_xx, -w_yy, -2*w_xy


def plot(size, index):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    X, Y = np.meshgrid(np.linspace(0, a, size), np.linspace(0, b, size))
    Z = np.zeros((size, size))

    for i in range(size):
        for j in range(size):
            STRAINS = np.array(strains(X[0, i], Y[j, 0], 10))
            NM = np.dot(ABD, STRAINS[:, None])
            Z[i,j] = NM[3+index]

    i,j = np.where(Z==Z.max())

    val = ['Mx', 'My', 'Mxy'][index]

    print('Max %s=%.4f at x,y=(%.4f,%.4f)' %(val, Z[i,j], Y[i,j], X[i,j]))



    ax.scatter(X, Y, Z)
    plt.show()
