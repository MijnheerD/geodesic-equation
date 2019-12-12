# Let's go!

from sympy import *
import numpy as np
from scipy.integrate import solve_ivp, odeint
from scipy.constants import G


class Spacetime:
    def __init__(self, g, symbols):
        self.metric = g
        self.__symbols = symbols
        self.diff = [self.diffFullMetric(mu) for mu in self.__symbols]

    def symbols(self): #Maybe use this to only give back the symbols used in the metric?
        return self.__symbols

    def diffFullMetric(self, mu):
        shape = self.metric.shape
        r = shape[0]
        c = shape[1]
        d = [x[:] for x in [[0] * r] * c]
        for row in range(r):
            for column in range(c):
                d[row][column] = diff(self.metric[row][column], mu)
        return d

    def diffMetric(self, pos, mu): #Mu given as number
        return diff(self.metric[pos[0], pos[1]], self.__symbols[mu])

    def chrisSymbol(self, up, down1, down2):
        return 0.5*sum([self.metric[up, m] * self.diffMetric([m, down1], down2)+self.metric[up, m] * self.diffMetric([m, down2], down1)-self.metric[up, m] * self.diffMetric([down1, down2], m) for m in range(len(self.__symbols))])

    def chrisSymbolFast(self, up, down1, down2):
        return 0.5*sum([self.metric[up][m] * self.diff[down2][m][down1] + self.metric[up][m] * self.diff[down1][m][down2] - self.metric[up][m] * self.diff[m][down1][down2] for m in range(len(self.__symbols))])


def geodesicEq(x, s, spacetime):
    if len(x) != 8:
        raise TypeError("Input vector must have length 8, not "+str(len(x)))
    if type(spacetime) != Spacetime:
        raise TypeError("Input spacetime must be of type Spacetime")
    u = [0]*8
    for i in range(4):
        u[i] = x[i+4]
    for i in range(4, 8):
        v = [spacetime.chrisSymbol(i-4, a, b)*u[a]*u[b] for a in range(4) for b in range(4)]
        p = -sum(v)
        u[i] = sp.lambdify(spacetime.symbols(), p, 'numpy')(x[0], x[1], x[2], x[3])
    return u


def solveGE(equation, xinit, ds, s0, s1, spacetime):
    s = np.arange(s0, s1+ds, ds)
    solution = solve_ivp(lambda s, x: equation(x, s, spacetime), [s0, s1], xinit, t_eval= np.arange(s0, s1+ds, ds), method='LSODA')
    return solution
  

def list_mul(lst, f):
    return [f * elem for elem in lst]


def list_sum(lst1, lst2):
    return [lst1[i] + lst2[i] for i in range(len(lst1))]


def rkStep(yn, f, tn, h, spacetime):
    k1 = [h * elem for elem in f(yn, tn, spacetime)]
    k2 = [h * elem for elem in f(list_sum(yn, list_mul(k1, 0.5)), tn + 0.5 * h, spacetime)]
    k3 = [h * elem for elem in f(list_sum(yn, list_mul(k1, 0.5)), tn + 0.5 * h, spacetime)]
    k4 = [h * elem for elem in f(list_sum(yn, k3), tn + h, spacetime)]
    return list_sum(yn, [(1.0 / 6.0) * elem for elem in (k1 + 2 * k2 + 2 * k3 + k4)])


def RK4(equation, xinit, ds, s0, s1, spacetime):
    print("Going into RK4")
    s = np.arange(s0, s1+ds, ds)
    sol = []
    x = xinit
    for t in s:
        print(t)
        x = rkStep(x, equation, t, ds, spacetime)
        sol.append(x[0:4])
    return sol
