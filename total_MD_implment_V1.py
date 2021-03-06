# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:15:32 2021

@author: chari
"""

#implement MD Code
import total_MD_code_V1 as F1
#import Lattice_Graphite as coor
import LAttice_SiO2 as coor
import numpy as np
import matplotlib.pyplot as plt
import time

def advance(pos, vel, mass, dt, disp, dist, rc, L,kvecs):
    """
    advance system according to velocity verlet

    args:
        pos (array): coordinates of particles
        val (array): velocities of particles
        mass (float): mass of particles
        dt (float): timestep by which to advance
        disp (array): displacement table
        dist (array): distance table
        rc (float): cutoff
        L (float): length of cubic box
    returns:
        array, array, array, array:
        new positions, new velocities, new displacement table,
        and new distance table
    """
    #print(F1.Ewald_force(pos,kvecs,alpha,V,L))
    accel = (F1.force(disp, dist, rc) + F1.Ewald_force(pos,kvecs,alpha,V,L)) / mass
    #move
    vel_half = vel + 0.5*dt*accel
    pos_new = pos #+ dt*vel_half
    pos_new[:,[0,1,2]] = pos[:,[0,1,2]] + dt*vel_half
    pos_new = F1.minimum_image(pos_new, L)
    #print(pos_new)
    disp_new = F1.displacement_table(pos_new, L)
    dist_new = np.linalg.norm(disp_new, axis=-1)
    #repeat force calculation for new pos
    accel = F1.force(disp_new, dist_new, rc) / mass
    #finish move
    vel_new = vel_half + 0.5*dt*accel
    return pos_new, vel_new, disp_new, dist_new

#initial Values
T = 1.0
L = 6.0
M = 1.0
cutoff = L/2
N = 44
V = L**3
alpha = (np.pi)*(N/(V**2))**(1/3)

steps = 5
timestep = 0.03

N_max = 4
Karr = F1.my_legal_kvecs(N_max,L)

num_bins = 14
dr_ = 0.25

#system
coordinates = coor.a
#Qs = np.random.uniform(-1,1,len(coordinates))
#Qs = np.array([Qs]).T
#coordinates = np.append(coordinates,Qs,axis=1)
velocities = F1.initial_velocities(N, M, T)

coorTEst = coordinates.copy()

#tables required to compute quantities like forces, energies
displacements = F1.displacement_table(coordinates[:,[0,1,2]], L)
distances = np.linalg.norm(displacements, axis=-1)

#advance and record energies
#it can also be useful to save coordinates and velocities
KE = []
PE = []
PE_EW = []
E = []
P = []
Temp = []
K_vecs = F1.my_legal_kvecs(N_max,L)
S_k = []

G_r = []
r_ = []

veloci = []
t0 = time.time()
for _ in range(steps):
    coordinates, velocities, displacements, distances = advance(coordinates,\
            velocities, M, timestep, displacements, distances, cutoff,\
            L,K_vecs)
    
    PE.append(F1.potential(distances, cutoff)+F1.Ewald_pot(coordinates,F1.my_legal_kvecs(N_max,L),alpha,V,L))
    PE_EW.append(F1.Ewald_pot(coordinates,F1.my_legal_kvecs(N_max,L),alpha,V,L))
    KE.append(F1.kinetic(M, velocities,N)[0])
    E.append(F1.potential(distances, cutoff) + F1.kinetic(M, velocities,N)[0])
    P.append(F1.kinetic(M, velocities,N)[1])
    Temp.append(F1.kinetic(M, velocities,N)[2])
    
    
    S_k.append(F1.my_calc_sk(F1.my_legal_kvecs(N_max,L),coordinates))
    
    veloci.append(velocities)
    
    G_r.append(F1.my_pair_correlation(F1.dist_ravel(distances),N,num_bins,dr_,L)[0])
    r_.append(F1.my_pair_correlation(F1.dist_ravel(distances),N,num_bins,dr_,L)[1])
    
t1 = time.time()
print(t1-t0)
























