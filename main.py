import time
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

    for i in range(50):
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

adam = pc.Particle([0,0,0],[0,0,-5e9])
sierra = pc.Particle([0,0,0],[0,0,5e9])

particles = [adam,sierra]

def init():
    for p in particle_array:
        ax.scatter(p.x_position,p.y_position,p.z_position, c='b', marker='o')

h = (1-0)/100
t = np.arange(0,1,h)
total = len(t)

time_array = []

start_time = time.time()
def animate(i):
#for instance in t:
    print(f"Iteration {i}/{total} starting...")
    iteration_time = time.time()
    mp.cla()
    ax.set_ylim(-box,box)
    ax.set_xlim(-box,box)
    ax.set_zlim(-box,box)
    temp = 0.0
    kinetic_energy = 0.0
    for p in particle_array:
        kinetic_energy += p.get_kinetic()
    Temprature = (2/3)*kinetic_energy/1.38e-23
    for p in particle_array:
        for q in particle_array:
            p.monte_carlo(q,Temprature,box)
        ax.scatter(p.x_position,p.y_position,p.z_position, c='b', marker='o')

    iteration_time_final = time.time()
    time_array.append(iteration_time_final - iteration_time)
    print(f"Iteration {i}/{total} Complete\n")
    print(f"Time elapsed:\t\t{time.time()-start_time}s")

ani = ap.FuncAnimation(fig,animate,t,init_func=init)
ani.save('mcAnimation.mp4', writer='ffmpeg', fps=60)

fig = mp.figure()
my_plot = fig.add_subplot()
mp.plot(time_array)
mp.show()
