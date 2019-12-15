from project import *
from scipy.constants import G
#import time

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


#Testing (same conditions as EinsteinPy)
s0 = 0
s1 = 5 #0.01
ds = 1 #0.3e-6
x = [0, 306.0, np.pi/2, -np.pi/6, 0, 0, 0, 1900]

#geodesicEq
#print(geodesicEq(x,2,ST2))

#solver
#u = rk4(geodesicEq, x, ds, s0, 2, ST2)
#print(u)

#u2 = solveGE(geodesicEq, x, ds, s0, s1, ST)
#print(u2.y)
#write_out([u2.t] + list(u2.y[1:4]), 'cart', 'solveGE_Minkowski.txt')

t3, y3 = rk4(geodesicEq, x, ds, s0, s1, ST)
x = []
y = []
z = []
for i in range(len(y3)):
    x.append(y3[i][1])
    y.append(y3[i][2])
    z.append(y3[i][3])
write_out([t3] + [x] + [y] + [z], 'cart', 'RK45_Minkowski.txt')
