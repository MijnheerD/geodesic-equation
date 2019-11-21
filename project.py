# Let's go!

from sympy import *
import numpy as np
from scipy.integrate import odeint


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
        return diff(self.metric[pos[0], pos[1]], self.__symbols[mu])

    def chrisSymbol(self, up, down1, down2):
        return 0.5*sum([self.metric[up, m] * self.diffMetric([m, down1], down2)+self.metric[up, m] * self.diffMetric([m, down2], down1)-self.metric[up, m] * self.diffMetric([down1, down2], m) for m in range(4)])

#I added an s argument in order to use odeint for solving the difeq
def geodesicEq(x, s, spacetime):
        if len(x)!=8:
            raise TypeError("Input vector must have length 8")
        if type(spacetime) != Spacetime:
            raise TypeError("Input spacetime must be of type Spacetime")
        u = [0]*8
        for i in range(4):
            u[i] = x[i+4]
        for i in range(4,8):
            v = [spacetime.chrisSymbol(i-4,a,b)*u[a]*u[b] for a in range(4) for b in range(4)]
            p = -sum(v)
            u[i] = p
        return u

#This function doesn't work yet as I don't know how odeint works with s as there is no s dependency in geodesicEq    
def solveGE(geodesicEq,xinit,ds,s0,s1,spacetime):
    s=np.arrange(s0,s1+ds,ds)
    solution=odeint(geodesicEq,xinit,s,args=(spacetime,))
    sol=[]
    for i in solution:
        sol.append([i[0:4]])
    return sol



#Examples of Metrics
t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])
ST = Spacetime(g)

g_schwartz=np.array([[-(1-1/x),0,0,0],[0,1/(1-1/x),0,0],[0,0,x**2,0],[0,0,0,x**2*np.sin(3)]])
ST2= Spacetime(g_schwartz)

#Test GeodesicEq
x = [1,1,2,3,3,4,3,5]
print(geodesicEq(x,ST2))
