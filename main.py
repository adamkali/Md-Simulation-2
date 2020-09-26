from sys import argv
from random import uniform
import matplotlib.pyplot as mp
import numpy as np
import particle_class as pc

def initialize(box):
    xyz = open("input.xyz", "w")

    for i in range(500):
        x, y, z = uniform(-box,box), uniform(-box,box), uniform(-box,box)
        print(x,y,z, file=xyz)

def init_particles(coords):

    # Finds the number of particles in the system
    size = coords.shape
    size = size[0]

    # Initialize a particle array that uses every single particle in the
    # box as an index.
    particle_array = []
    for p in range(size):
        coord = np.array(coords[p][:])
        particle_array.append(pc.Particle(coord))

    return particle_array

def init_eos(particles,U,box):
    kinetic_energy = 0.0
    for p in particles:
        kinetic_energy += p.get_kinetic()
    Temprature = (2/3)*kinetic_energy/1.38e-23

    potential = 0.0
    for p in particles:
        for q in particles:
            potential += p.leonard_jones(q)
    U.append(potential)

    N = len(particles)
    kb,a,b = 1.3806e-23,553.6,0.0305
    Temprature_prime = Temprature * 358.38
    Avo = 6.022e23
    pv_ideal = N*kb*Temprature_prime
    n = N/Avo
    V = (2*box*1e-10)**3
    pv_ideal = pv_ideal/(V - n*b)
    return pv_ideal - (a*(n/V)**2), Temprature


def eos(particles,Temp,box):

    V = (2*box*1e-10)**3
    n = len(particles)/6.022e23
    R, a, b =  8.3145, 553.6, 0.0305
    Temprature_prime = Temp * 358.38
    pv_ideal = n*R*Temprature_prime
    pv_ideal = pv_ideal/(V - n*b)

    return pv_ideal - (a*(n/V)**2)

def main():

    box = float(argv[1])

    initialize(box)

    h = (1-0)/100
    t = np.arange(0,1,h)

    U = []
    pressure = []

    # Load in input.xy to initalize the system
    coordinates = np.loadtxt("input.xyz")
    particle_array = init_particles(coordinates)
    P, Temp = init_eos(particle_array,U,box)
    pressure.append(P)

    for instance in range(len(t)):
        potential = 0.0
        for p in particle_array:
            for q in particle_array:
                p.monte_carlo(q,Temp,box)
                potential += p.leonard_jones(q)
        U.append(potential)
        pressure.append(eos(particle_array,Temp,box))




    print(f"Average Temprature\t\t\t{np.mean(Temp)}\nAverage LJ Energy\t\t\t{np.mean(U)}")

    mp.plot(pressure)
    mp.plot(U)
    mp.show()

if __name__ == "__main__":
    main()
