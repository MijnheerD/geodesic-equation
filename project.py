# Let's go!

from sympy import *
import numpy as np
from scipy.integrate import solve_ivp, odeint


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
            raise TypeError("Input vector must have length 8, not "+str(len(x)))
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

#This function gives A type error for line solution=odeint(geodesicEq,xinit,s,args=(spacetime,)) when calling it   
def solveGE(geodesicEq,xinit,ds,s0,s1,spacetime):
    s=np.arange(s0,s1+ds,ds)
    solution=odeint(geodesicEq,xinit,s,args=(spacetime,))
    sol=[]
    for i in solution:
        sol.append([i[0:4]])
    return sol
#Let's try with RK4=> still an error but now it says "
#TypeError: can't multiply sequence by non-int of type 'float'
#and the problem lays in k1 = h * f(yn,tn,spacetime)

def list_mul(lst, f):
    return [f * elem for elem in lst]

def list_sum(lst1, lst2):
    return [lst1[i] + lst2[i] for i in range(len(lst1))]

def rkStep(yn, f, tn, h,spacetime):
        k1 = [h * elem for elem in f(yn,tn,spacetime)]
        k2 = [h * elem for elem in f(list_sum(yn, list_mul(k1, 0.5)), tn + 0.5 * h, spacetime)]
        k3 = [h * elem for elem in f(list_sum(yn, list_mul(k1, 0.5)), tn + 0.5 * h, spacetime)]
        k4 = [h * elem for elem in f(list_sum(yn, k3), tn + h, spacetime)]
        return list_sum(yn, [(1.0 / 6.0) * elem for elem in (k1 + 2 * k2 + 2 * k3 + k4)])

def RK4(geodesicEq,xinit,ds,s0,s1,spacetime,h):
    s=np.arange(s0,s1+ds,ds)
    sol=[]
    x = xinit
    for t in s:
        x = rkStep(x,geodesicEq,t,h,spacetime)
        sol.append(x[0:4])
    return sol



#Examples of Metrics
t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

print(0.5*t)

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])
ST = Spacetime(g)

g_schwartz=np.array([[-(1-1/x),0,0,0],[0,1/(1-1/x),0,0],[0,0,x**2,0],[0,0,0,x**2*np.sin(3)]])
ST2= Spacetime(g_schwartz)


#Testing
s0 = 1
s1 = 5
ds = 1
s = np.arange(s0,s1+ds,ds)
x = [1,1,2,3,3,4,3,5]

#geodesicEq
#print(geodesicEq(x,2,ST2))

#solver
u=RK4(geodesicEq,x,ds,s0,s1,ST2,0.5)
print(u)

#sol = []
#solution = solve_ivp(lambda t, y: geodesicEq(y,t,ST), (s0,s1), x,args=(ST,))
#for i in solution:
#    sol.append([i[0:4]])
#print(sol)