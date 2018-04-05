import numpy as np
from sympy import *


x,y = var('x y')
a,b = var('a b')
m,n = var('m n')

q = (4500*x*(a-x)**2)/(a*b)

equ = q * sin(m*pi*x/a) * sin(n*pi*y/b)

subs = {
    a: 1,
    b: 1,
    m: 2,
    n: 3,
    }

Q = integrate( integrate(equ, (x,0,a)), (y,0,b))
