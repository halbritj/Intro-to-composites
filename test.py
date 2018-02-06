import numpy as np

a = np.random.rand(9).reshape(3,3)
b = np.random.rand(9).reshape(3,3)

A = np.linalg.inv(np.dot(a,b))
B = np.dot(np.linalg.inv(a), np.linalg.inv(b))
