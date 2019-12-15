#from project import read_out
from visual import *
#from scipy.constants import G
from math import sin, cos
G = 6.67 * 10**(-11)
c = 2.998 * 10**8


def read_out(filename):
    f = open(filename, "r")
    line = f.readline()
    coord = line.split()[0]
    line = f.readline()
    t = [float(elem) for elem in line.split()]
    line = f.readline()
    r = [float(elem) for elem in line.split()]
    line = f.readline()
    theta = [float(elem) for elem in line.split()]
    line = f.readline()
    phi = [float(elem) for elem in line.split()]
    f.close()

    if coord == 'spherical':
        x = []
        y = []
        z = []
        for i in range(len(theta)):
            y.append(r[i] * sin(theta[i]) * sin(phi[i]))
            x.append(r[i] * sin(theta[i]) * cos(phi[i]))
            z.append(r[i] * cos(theta[i]))
        return (t, x, y, z)
    elif coord == 'cart':
        return (t, r, theta, phi)
    else:
        raise NameError("Coordinate system unknown")



def animate_schwarzschild_geodesic(file, M, rt):
    t, x, y, z = read_out(file)

    A = vector(x[0], y[0], z[0])
    B = vector(x[1], y[1], z[1])
    C = cross(A,B)

    scene = display(center = (0,0,0), forward = -norm(C), fullscreen = True, exit = False, background = (0.3,0.3,0.3))
    
    r_s = 2*G*M/(c**2)
    f = frame()
    particle = sphere(frame = f, pos = vector(x[0], y[0], z[0]), radius = 5, color = color.green, make_trail = True)
    black_hole = sphere(frame = f, pos = vector(0,0,0), radius = 10, color = color.black)
    
    for i in range(len(t)):
        rate(rt)
        particle.pos = vector(x[i], y[i], z[i])

    particle.color = color.red


def animate_minkowski_geodesic(file, rt):
    t, x, y, z = read_out(file)

    scene = display(center = vector(x[0], y[0], z[0]), forward = (-1, 0, 0), range = (5,5,5), fullscreen=True, exit=False, background=(0.3,0.3,0.3))

    particle = sphere(canvas = scene, pos=vector(x[0], y[0], z[0]), radius = 0.1, color=color.green, make_trail=True)

    for i in range(len(t)):
        rate(rt)
        particle.pos = vector(x[i], y[i], z[i])

    particle.color = color.red


#animate_minkowski_geodesic("solveGE_Minkowski.txt", 1000)
#animate_minkowski_geodesic("RK4_Minkowski.txt", 2000)
#animate_schwarzschild_geodesic("solveGE_phi.txt", 5.972e24, 10000)
animate_schwarzschild_geodesic("solveGE_theta.txt", 5.972e24, 900)
#animate_schwarzschild_geodesic("solveEinstein.txt", 5.972e24, 100)
#animate_schwarzschild_geodesic("RK4.txt", 5.972e24, 5000)
