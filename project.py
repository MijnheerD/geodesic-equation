# Let's go!

from sympy import *
import numpy as np


#Should make a class "Spacetime" with property g

class Spacetime:
    def __init__(self, g):
        self.metric = g
        self.__symbols = [t, x, y, z]

    def diffFullMetric(self, mu):  # Overkill for Christoffel symbols, might be useful later on
        shape = self.metric.shape
        d = np.full(shape, 0)
        r = len(g)
        c = len(g[0])
        for row in range(r):
            for column in range(c):
                d[row, column] = diff(g[row, column], self.__symbols[mu])
        return d

    def diffMetric(self, pos, mu): #Mu given as number
        return diff(g[pos[0], pos[1]], self.__symbols[mu])

    def chrisSymbol(self, up, down1, down2):
        return 0.5*sum([g[up, m] * self.diffMetric([m, down1], down2)+g[up, m] * self.diffMetric([m, down2], down1)-g[up, m] * self.diffMetric([down1, down2], m) for m in range(4)])


t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])

ST = Spacetime(g)
print(ST.chrisSymbol(1,1,1))