import random as rand
from matplotlib.pyplot import plot
from math import sqrt
import numpy as np

class Particle():

    #initiate the class
    def __init__(self, arr, vel_arr=None):
        self.x_position = arr[0]
        self.y_position = arr[1]
        self.z_position = arr[2]
        if vel_arr == None:
            self.x_velocity = rand.uniform(-4.04,4.04)
            self.y_velocity = rand.uniform(-4.04,4.04)
            self.z_velocity = rand.uniform(-4.04,4.04)
        else:
            self.x_velocity = vel_arr[0]
            self.y_velocity = vel_arr[1]
            self.z_velocity = vel_arr[2]

        # Define interaction constants for later use
        self.epsilon = 4.91151e-21
        self.sigma = 2.725


    # Define __str__ for debugging
    def __str__(self):
        return f"Particle:\nx:\t({self.x_position},{self.x_velocity})\ny:\t({self.y_position},{self.y_velocity})\nz:\t({self.z_position},{self.z_velocity})\n"

    # Define equate to tell if a particle is itself
    def __eq__(self, other):
        return self.x_position == other.x_position \
                and self.x_velocity == other.x_velocity \
                and self.y_position == other.y_position \
                and self.y_velocity == other.y_velocity \
                and self.z_position == other.z_position \
                and self.z_velocity == other.z_velocity

    def get_poss(self):
        return self.x_position,self.y_position,self.z_position

    def get_vels(self):
        return self.x_velocity,self.y_velocity,self.z_velocity

    def leonard_jones(self, other):

        # Account for when there is a particle interating with itself
        if self == other:
            return 0

        # Now for the meat of the model
        else:
            X = (self.x_position - other.x_position)
            Y = (self.y_position - other.y_position)
            Z = (self.z_position - other.z_position)
            R = sqrt(X**2 + Y**2 + Z**2)
            return 4*self.epsilon*((self.sigma/R)**12 - (self.sigma/R)**6)

    def get_kinetic(self):
        mass = 30.113e-27 # kg
        return 0.5*mass*(self.x_velocity**2+self.y_velocity**2+self.z_velocity**2)

    def monte_carlo(self, other, Temp, box):
        random_number= rand.uniform(0,1)
        kb = 1.3806e-23
        random_array = [rand.uniform(-box,box),
                rand.uniform(-box,box),
                rand.uniform(-box,box)]
        ghost = Particle(random_array,
                [self.x_velocity,self.y_velocity,self.z_velocity])
        delta_pot = ghost.leonard_jones(other) - self.leonard_jones(other)
        if delta_pot <= 0:
            self.x_position = ghost.x_position
            self.y_position = ghost.y_position
            self.z_position = ghost.z_position
        elif np.exp(-delta_pot/(kb*Temp)) < random_number:
            self.x_position = ghost.x_position
            self.y_position = ghost.y_position
            self.z_position = ghost.z_position
