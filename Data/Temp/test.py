#!/usr/bin/env python

import copy as cp
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker as tk
from scipy import ndimage as ni
from scipy import interpolate as it

# (1) - Interpolation test.
'''
v = 2
m = 1.7

a = np.arange(0, 12.001, 1)
b = np.sin(v*a)
rep = it.splrep(a, b, s=0)

newa = np.arange(0, 12.001, 0.17)
newb = it.splev(newa, rep, der=0)

cx = np.arange(0.3, 12.001, m)
cy = it.splev(cx, rep, der=0)

f = np.arange(0, 12.001, 0.01)

tx = np.arange(0.3, 12.001, m)
ty = []
for i in range(len(tx)):
    ty.append([])
for i in range(len(newa)):
    ty[np.argmin(np.absolute(newa[i]-tx))].append(newb[i])
for i in range(len(ty)):
    ty[i] = np.mean(ty[i])

pl.plot(a, b, 'bx', newa, newb, 'co', f, np.sin(v*f), 'r:', tx, ty, 'kd', cx, cy, 'gd')
pl.grid(True)
pl.show()
pl.close()
'''

# (2) - Array test.
'''
a = np.array([[1, 2], [7, 3, 4, 5], [-1], [0, 0, 1]])

print(a)

a = np.insert(a, 1, [7, 8])

print(a)

b = np.array([2, 4, 5, 6, 8, 9, 13, 17])

c = np.array([3, 4, 6])

print(b[c])
'''

# (3) - Copy test.
'''
a = [7.0, np.array(['3', 3, 3.0]), 'a', 3.5]

b = cp.deepcopy(a)

c = a

d = a

print(type(a), type(b), type(c), type(d), '\n', a, b, c, d)

c[1][1] = 5

print(type(a), type(b), type(c), type(d), '\n', a, b, c, d)

e = ['e', 4, 7, 8.5]

f = np.copy(e)

g = cp.deepcopy(e)

print(type(e), type(f), type(g))
'''

# (4) Numpy test.
'''
a = [[1, 2, 3], [4, 5], [6, -2, -3, 0, 1]]

b = np.array(a)

c = cp.deepcopy(list(b))

print(type(a), type(a[0]), type(b), type(b[0]), type(c), type(c[0]))

m = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]

print(m)

print(np.array(m))

print(np.insert(m, len(m), 0, axis=0))
'''

# (5) Location test.
'''
a = [0, 2, 2.0, 4/2, 5, 7, 2.0, -2, 2, 2.0/1, 3, -8, -1.2]

b = np.delete(range(len(a)), np.nonzero(np.array(a)-2))

print(a, '\n', b)
'''

# (6) Convolution test.
'''
mname = '/home/diogo/Documents/Thesis/Data/Models/spec_1-1-1-1-1-1.dat'
model = np.flip(np.swapaxes(np.genfromtxt(mname, dtype=float), 0, 1), 1)
conv = ni.filters.gaussian_filter(model[1], 3/(2*(2*np.log(2))**0.5))

pl.plot(model[0], model[1], 'r:', model[0], conv, 'k-.')
pl.grid(True)
pl.show()
pl.close()
'''
