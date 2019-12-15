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
        y = [r[i] * sin(theta[i]) * sin(phi[i]) for i in range(len(theta))]
        x = [r[i] * sin(theta[i]) * cos(phi[i]) for i in range(len(theta))]
        z = [r[i] * cos(theta[i]) for i in range(len(theta))]
        return (t, x, y, z)
    elif coord == 'cart':
        return (t, r, theta, phi)
    else:
        raise NameError("Coordinate system unknown")


def animate_geodesic(file, M, rt):
    t, x, y, z = read_out(file)

    scene = display(fullscreen = True, exit = False, background = color.white)
    
    r_s = 2*G*M/(c**2)
    f = frame()
    particle = sphere(frame = f, pos = vector(x[0], y[0], z[0]), radius = 5, color = color.green, make_trail = True)
    black_hole = sphere(frame = f, pos = vector(0,0,0), radius = 1, color = color.black)
    
    for i in range(len(t)):
        rate(rt)
        particle.pos = vector(x[i], y[i], z[i])

    particle.color = color.red


animate_geodesic("solveGE.txt", 5.972e24, 5000)
animate_geodesic("solveEinstein.txt", 5.972e24, 100)
