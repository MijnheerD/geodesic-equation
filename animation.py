from project import *
from visual import *
#from scipy.constants import G
from math import sin, cos
G = 6.67 * 10**(-11)


def animate_geodesic(file, M):
    t, x, y, z = read_out(file)
    
    r_s = 2*G*M
    particle = sphere(pos = vector(x[0], y[0], z[0]), radius = 1, color = color.green)
    black_hole = sphere(pos = vector(0,0,0), radius = r_s, color = color.red)
    
    for i in range(len(t)):
        rate(1)
        particle.pos = vector(x[i], y[i], z[i])
