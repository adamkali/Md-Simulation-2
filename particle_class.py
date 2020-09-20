import random as rand
from matplotlib.pyplot import plot
from math import sqrt
import numpy as np

class Particle():


        #Define mass to determine Average Kinetic Energy and Temprature
    mass = 2.9915e-20 # kg
    c = 3e8

    #initiate the class
    def __init__(self, arr, vel_arr=None):
        self.x_position = arr[0]
        self.y_position = arr[1]
        self.z_position = arr[2]
        if vel_arr == None:
            self.x_velocity = rand.uniform(-10,10)
            self.y_velocity = rand.uniform(-10,10)
            self.z_velocity = rand.uniform(-10,10)
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

    def leonard_jones(self, other, h1=0, h2=0, h3=0):

        # Account for when there is a particle interating with itself
        if self == other:
            return 0

        # Now for the meat of the model
        else:
            X = (self.x_position - other.x_position) + h1
            Y = (self.y_position - other.y_position) + h2
            Z = (self.z_position - other.z_position) + h3
            R = sqrt(X**2 + Y**2 + Z**2)
            return 4*self.epsilon*((self.sigma/R)**12 - (self.sigma/R)**6)


    def force(self,other):


        hx, hy, hz = ((self.x_position-other.x_position + .001)/100),\
                ((self.y_position-other.y_position + 0.001)/100 ), \
                ((self.z_position-other.z_position + .001)/100)

        F_x = (self.leonard_jones(other,h1=hx)-self.leonard_jones(other,h1=-hx))/(2*hx)
        F_y = (self.leonard_jones(other,h2=hy)-self.leonard_jones(other,h2=-hy))/(2*hy)
        F_z = (self.leonard_jones(other,h3=hy)-self.leonard_jones(other,h3=-hz))/(2*hz)
        return np.array([F_x,F_y,F_z], float)

    def verlet(self, other,step_size):


        if self.get_poss() == other.get_poss():

            print("-"*40 + "collision" + "-"*40)
            self.x_position += -self.x_velocity*step_size
            self.y_position += -self.y_velocity*step_size
            self.z_position += -self.z_velocity*step_size
            other.x_position += -other.x_velocity*step_size
            other.y_position += -other.y_velocity*step_size
            other.z_position += -other.z_velocity*step_size

        verlet_force = self.force(other)
        half_step = 0.5*step_size*verlet_force
        vx_half = self.x_velocity + half_step[0]
        vy_half = self.y_velocity + half_step[1]
        vz_half = self.z_velocity + half_step[2]
        self.x_position += step_size*vx_half
        self.y_position += step_size*vy_half
        self.z_position += step_size*vz_half
        ghost = Particle(np.array([self.x_position,self.y_position,self.z_position], float),
                [vx_half, vx_half, vz_half])
        k = ghost.force(other)
        self.x_velocity += step_size*k[0]
        self.y_velocity += step_size*k[1]
        self.z_velocity += step_size*k[2]
        if self.x_velocity >= 50:
            print("Problem is here...")
            self.x_velocity = 50
        if self.y_velocity >= 50:
            print("Problem is here...")
            self.y_velocity = 50
        if self.z_velocity >= 50:
            print("Problem is here...")
            self.z_velocity = 50

    def boundaries(self,box):
        if self.x_position >= box:
            if self.y_position >= box:
                if self.z_position >= box:
                    return box,box,box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return box,box,-box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                else:
                    return box,box,self.z_position,-self.x_velocity,-self.y_velocity,self.z_velocity
            elif self.y_position <= -box:
                if self.z_position >= box:
                    return box,-box,box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return box,-box,-box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                else:
                    return box,-box,self.z_position,-self.x_velocity,-self.y_velocity,self.z_velocity
            else:
                if self.z_position >= box:
                    return box,self.y_position,box,-self.x_velocity,self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return box,self.y_position,-box,-self.x_velocity,self.y_velocity,-self.z_velocity
                else:
                    return box,self.y_position,self.z_position,-self.x_velocity,self.y_velocity,self.z_velocity

        elif self.x_position <= -box:
            if self.y_position >= box:
                if self.z_position >= box:
                    return -box,box,box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return -box,box,-box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                else:
                    return -box,box,self.z_position,-self.x_velocity,-self.y_velocity,self.z_velocity
            elif self.y_position <= -box:
                if self.z_position >= box:
                    return -box,-box,box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return -box,-box,-box,-self.x_velocity,-self.y_velocity,-self.z_velocity
                else:
                    return -box,-box,self.z_position,-self.x_velocity,-self.y_velocity,self.z_velocity
            else:
                if self.z_position >= box:
                    return -box,self.y_position,box,-self.x_velocity,self.y_velocity,-self.z_velocity
                elif self.z_position <= -box:
                    return -box,self.y_position,-box,-self.x_velocity,self.y_velocity,-self.z_velocity
                else:
                    return -box,self.y_position,self.z_position,-self.x_velocity,self.y_velocity,-self.z_velocity

        elif self.y_position >= box:
            if self.z_position >= box:
                return self.x_position,box,box,self.x_velocity,-self.y_velocity,-self.z_velocity
            elif self.z_position <= -box:
                return self.x_position,box,-box,self.x_velocity,-self.y_velocity,-self.z_velocity
            else:
                return self.x_position,box,self.z_position,self.x_velocity,-self.y_velocity,self.z_velocity

        elif self.y_position <= -box:
            if self.z_position >= box:
                return self.x_position,-box,box,self.x_velocity,-self.y_velocity,-self.z_velocity
            elif self.z_position <= -box:
                return self.x_position,-box,-box,self.x_velocity,-self.y_velocity,-self.z_velocity
            else:
                return self.x_position,-box,self.z_position,self.x_velocity,-self.y_velocity,self.z_velocity

        elif self.z_position >= box:
            return self.x_position,self.y_position,box,self.x_velocity,self.y_velocity,-self.z_velocity

        elif self.z_position <= -box:
            return self.x_position,self.y_position,-box,self.x_velocity,self.y_velocity,-self.z_velocity

        else:
            return self.x_position,self.y_position,self.z_position,self.x_velocity,self.y_velocity,self.z_velocity



