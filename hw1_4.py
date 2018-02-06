import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


def L(plane):
    if plane == '12': return sp.diag(1,1,-1)
    if plane == '23': return sp.diag(-1,1,1)
    if plane == '31': return sp.diag(1,-1,1)
    



s1, s2, s3, s4, s5, s6 = sp.var('s1 s2 s3 s4 s5 s6')

Stress = sp.Matrix([
    [s1, s6, s5],
    [s6, s2, s4],
    [s5, s4, s3]])
