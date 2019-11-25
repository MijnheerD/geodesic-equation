from project import *

#Examples of Metrics
t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
theta=Symbol('theta')
r=Symbol('r')
phi=Symbol('phi')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])
ST = Spacetime(g)

g_schwartz=np.array([[-(1-1/r),0,0,0],[0,1/(1-1/r),0,0],[0,0,r**2,0],[0,0,0,x**2*sin(theta)]]
ST2= Spacetime(g_schwartz)
#alternative definition: rs=1
#g_schwartz=np.array([[-(1-1/r),0,0,0],[0,1/(1-1/r),0,0],[0,0,r**2,0],[0,0,0,x**2*sin(theta)]])
#ST2= Spacetime(g_schwartz,t,r,theta,phi)
#gives christoffelsymbol 0,0,0,...                  


#Testing
s0 = 1
s1 = 5
ds = 1
x = [2,3,4,5,2,3,4,5]

#geodesicEq
#print(geodesicEq(x,2,ST2))

#solver
u = RK4(geodesicEq,x,ds,s0,s1,ST2)
print(u)

u2 = solveGE(geodesicEq,x,ds,s0,s1,ST2)
print(u2)
