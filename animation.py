from project import *
from vpython impotrt *
from scipy.constants import G


t = Symbol('t')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

g = np.array([[t,0,0,0],[0,x,0,0],[0,0,y,0],[0,0,0,z]])
ST = Spacetime(g)

g_schwartz=np.array([[-(1-1/x),0,0,0],[0,1/(1-1/x),0,0],[0,0,x**2,0],[0,0,0,x**2*np.sin(3)]])
ST2= Spacetime(g_schwartz)

s0 = 1
s1 = 5
ds = 1
u0 = [2,3,4,5,2,3,4,5]

u = np.array(RK4(geodesicEq,u0,ds,s0,s1,ST2))
u2 = np.array(solveGE(geodesicEq,u0,ds,s0,s1,ST2))


def animate_geodesic(file, M):
    
    #read data from file
    
    t = u2.t
    r = u2.y[1]
    theta = u2.y[2]
    phi = u2.y[3]
    
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    r_s = 2*G*M
    particle = sphere(pos = vector(x[0], y[0], z[0]), radius = 1, color = color.green)
    black_hole = sphere(pos = vector(0,0,0), radius = r_s, color = color.red)
    
    for i in range(len(x)):
        rate(10)
        particle.pos = vector(x[i], y[i], z[i])
