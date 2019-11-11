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

#Should make a class "Spacetime" with property g

#No need for diffMetric?

class Spacetime:
    def __init__(self, g):
        self.metric = g
        self.symbols = [t, x, y, z]

    def chrisSymbol(self, up, down1, down2):
        t1 = sum([g[up, m] * diff(g[m, down1], down2) for m in self.symbols])
        t2 = sum([g[up, m] * diff(g[m, down2], down1) for m in self.symbols])
        t3 = sum([g[up, m] * diff(g[down1, down2], m) for m in self.symbols])
        return 0.5*(t1+t2-t3)

ST = Spacetime(g)
print(ST.chrisSymbol(0,1,2))