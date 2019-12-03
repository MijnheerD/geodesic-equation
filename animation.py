from project import *
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


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

data = [u,u2]

xl = [u[i][1] for i in range(5)]
yl = [u[i][2] for i in range(5)]
zl = [u[i][3] for i in range(5)]

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0:1,1], dat[0:1,2], dat[0:1,3])[0] for dat in data]
#ax.plot(xl, yl, zl)

# Setting the axes properties
ax.set_xlim3d([0.0, 1000.0])
ax.set_xlabel('X')

ax.set_ylim3d([0.0, 1000.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 1000.0])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
#line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),
#                                   interval=50, blit=False)

plt.show()
