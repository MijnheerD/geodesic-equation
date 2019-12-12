from project import *
import time

#Examples of Metrics
t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

theta = Symbol('theta')
r = Symbol('r')
phi = Symbol('phi')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])
ST = Spacetime(g, [t,x,y,z])

def schwartzmetric(M):
    return np.array([[-(1-2*G*M/r*2*10**8),0,0,0],[0,1/(1-2*G*M/r),0,0],[0,0,r**2,0],[0,0,0,r**2*sin(theta)]])

g_schwartz = schwartzmetric(5.972e24)
ST2 = Spacetime(g_schwartz, [t, r,theta,phi])

#This shows a factor of 10 improvement in runtime
'''
t0 = time.time_ns()
[ST2.chrisSymbol(0, a, b) for a in range(4) for b in range(4)]
t1 = time.time_ns()
[ST2.chrisSymbolFast(0, a, b) for a in range(4) for b in range(4)]
t2 = time.time_ns()
print(t1-t0)
print(t2-t1)
'''

#Testing (same conditions as EinsteinPy)
s0 = 0
s1 = 0.01
ds = 0.3e-6
x = [0, 306.0, np.pi/2, -np.pi/6, 0, 0, 0, 1900]

#geodesicEq
#print(geodesicEq(x,2,ST2))

#solver
#u = RK4(geodesicEq, x, ds, s0, 2, ST2)
#print(u)

#t0 = time.time_ns()
#u2 = solveGE(geodesicEq, x, ds, s0, s1, ST2)
#t1 = time.time_ns()
#print(u2)
#print("SolveGE took "+str(t1-t0)+" ns")

t0 = time.time_ns()
t3, y3 = rk45(geodesicEq, x, ds, s0, s1, ST2)
t1 = time.time_ns()
print(y3)
print("rk45 took "+str(t1-t0)+" ns")
