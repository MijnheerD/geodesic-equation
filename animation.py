#from project import *
from visual import *
#from scipy.constants import G
from math import sin, cos
G = 6.67 * 10**(-11)

def write_out(data, coord, filename):
    f = open(filename, "w")
    #Which coordinate system are we using for the data?
    f.write(coord)
    f.write("\n")
    #Assume data structure
    t = data[0]
    x = data[1]
    y = data[2]
    z = data[3]
    for item in t:
        f.write("%f " % item)
    f.write("\n")
    for item in x:
        f.write("%f " % item)
    f.write("\n")
    for item in y:
        f.write("%f " % item)
    f.write("\n")
    for item in z:
        f.write("%f " % item)


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

def animate_geodesic(file, M):
    t, x, y, z = read_out(file)
    
    r_s = 2*G*M
    particle = sphere(pos = vector(x[0], y[0], z[0]), radius = 1, color = color.green)
    black_hole = sphere(pos = vector(0,0,0), radius = r_s, color = color.red)
    
    for i in range(len(t)):
        rate(1)
        particle.pos = vector(x[i], y[i], z[i])

write_out([[0,1,2,3,4],[1,2,3,4,5],[1,2,3,4,5],[1,1,1,1,1]], "cart", "test.txt")
animate_geodesic("test.txt", 10**11)
