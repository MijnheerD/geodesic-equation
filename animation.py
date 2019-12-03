from project import *
from vpython import *
from scipy.constants import G


def write_out(data, filename):
    f = open(filename, "w")
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


def animate_geodesic(file, M):
    f = open(file, "r")
    line = f.readline()
    t = [float(elem) for elem in line.split()]
    line = f.readline()
    r = [float(elem) for elem in line.split()]
    line = f.readline()
    theta = [float(elem) for elem in line.split()]
    line = f.readline()
    phi = [float(elem) for elem in line.split()]
    f.close()

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    r_s = 2*G*M
    particle = sphere(pos = vector(x[0], y[0], z[0]), radius = 1, color = color.green)
    black_hole = sphere(pos = vector(0,0,0), radius = 1, color = color.red)
    
    for i in range(len(t)):
        rate(1)
        particle.pos = vector(x[i], y[i], z[i])

write_out([[0,1,2,3,4],[1,2,3,4,5],[1,2,3,4,5],[1,1,1,1,1]], "test.txt")
animate_geodesic("test.txt", 10**11)
