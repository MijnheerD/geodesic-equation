# Let's go!

from sympy import *
import numpy as np

t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])

def diffMetric(g,mu):
    shape = g.shape
    d = np.full(shape, 0)
    r = len(g)
    c = len(g[0])
    for row in range(r):
        for column in range(c):
            d[row, column] = diff(g[row,column],mu)
    return d

print(diffMetric(g,y))