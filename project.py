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
#test on Schwarzschild metric
#put rs=1, r=x theta=3
g_schwartz=np.array([[-(1-1/x),0,0,0],[0,1/(1-1/x),0,0],[0,0,x**2,0],[0,0,0,x**2*np.sin(3)]])

print(g[1,1])
ST = Spacetime(g)
ST2= Spacetime(g_schwartz)
print(ST.chrisSymbol(1,1,1))
print(ST2.chrisSymbol(3,3,3))
#Warning the last print gives us 0.5z, this doesn't seem to be normal
