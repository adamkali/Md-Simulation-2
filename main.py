from sys import argv
from random import uniform
import matplotlib.pyplot as mp
import numpy as np
import particle_class as pc
import matplotlib.animation as ap
from mpl_toolkits.mplot3d import Axes3D

"""====================================================================
Initialization
===================================================================="""

#initialize the box
box = float(argv[1])

def initialize(box):
    xyz = open("input.xyz", "w")

    for i in range(100):
        x, y, z = uniform(-box,box), uniform(-box,box), uniform(-box,box)
        print(x,y,z, file=xyz)

initialize(box)

# Load in input.xy to initalize the
coordinates = np.loadtxt("input.xyz")

# Finds the number of particles in the system
size = coordinates.shape
size = size[0]

# Initialize a particle array that uses every single particle in the
# box as an index.
particle_array = []
for p in range(size):
    coord = np.array(coordinates[p][:])
    particle_array.append(pc.Particle(coord))

fig = mp.figure()
ax = fig.add_subplot(111, projection='3d')

adam = pc.Particle([0,0,0],[0,0,-10])
sierra = pc.Particle([0,0,0],[0,0,10])

particles = [adam,sierra]

def init():
    for p in particle_array:
        ax.scatter(p.x_position,p.y_position,p.z_position, c='b', marker='o')

h = (1-0)/100
t = np.arange(0,10,h)

def animate(i):

    mp.cla()
    ax.set_ylim(-box,box)
    ax.set_xlim(-box,box)
    ax.set_zlim(-box,box)

    for p in particle_array:
        for i in particle_array:
            if p != i:
                p.verlet(i,h)

        p.x_position, p.y_position, p.z_position, p.x_velocity, p.y_velocity, p.z_velocity = p.boundaries(box)
        ax.scatter(p.x_position,p.y_position,p.z_position, c='b', marker='o')



ani = ap.FuncAnimation(fig,animate,t,init_func=init)
ani.save('myAnimation.mp4', writer='ffmpeg', fps=60)



