# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 17:14:32 2021

@author: chari
"""

##NEw Implement

#implement MD Code
import total_MD_code_V2 as F1
#import Lattice_Graphite as coor
import LAttice_SiO2 as coor
import numpy as np
import matplotlib.pyplot as plt
import time

##Initial Values
#Coordinates
initial_pos = coor.a 

L = coor.dim
T = 300
cutoff = L/2
N = np.shape(initial_pos)[0]
V = L**3
alpha = (np.pi)*(N/(V**2))**(1/3)
#
steps = 10
timestep =1e-12 
#
N_max = 4
Karr = F1.my_legal_kvecs(N_max,L)
#
num_bins = 14
dr_ = 0.25

NC = 9/7/10**9
Efield = np.array([NC,0,0])


#########
#Advance function
def advance(pos, vel, dt, disp, dist, rc, L,kvecs,E):
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
    accel = (F1.forceLJ(pos,rc,L) + F1.Ewald_force(pos,kvecs,alpha,V,L) + F1.external_force(pos,E))
    #print(accel)
    for i in range(len(accel)):
        m = pos[i][3]*(1/1000)
        accel[i] = accel[i]/m
    #print(accel)
    #move
    vel_half = vel + 0.5*dt*accel
    for i in range(0,np.shape(vel_half)[0]):
        if pos[i][4] == 1:
            vel_half[i][0:3] = vel[i][0:3] + 0.5*dt*accel[i][0:3]
        else:
            vel_half[i][0:3] = 0.0
    pos_new = F1.minimum_image(pos,L) #+ dt*vel_half
    for x in range(0,pos_new.shape[0]):
        if pos_new[x][4] == 1:  
            pos_new[x][0:3]= F1.minimum_image(pos[x,[0,1,2]] + dt*vel_half[x],L)
    #pos_new[:,[0,1,2]] = pos[:,[0,1,2]] + dt*vel_half
    #pos_new = F1.minimum_image(pos_new, L)
    #print(pos_new)
    disp_new = F1.displacement_table(pos_new, L)
    dist_new = np.linalg.norm(disp_new, axis=-1)
    #repeat force calculation for new pos
    accel = (F1.forceLJ(pos_new,rc,L) + F1.Ewald_force(pos_new,kvecs,alpha,V,L)+ F1.external_force(pos,E))
    for i in range(len(accel)):
        m = pos[i][3]*(1/1000)
        accel[i] = accel[i]/m
    #print(accel)
    #finish move
    vel_new = vel_half + 0.5*dt*accel
    #print(vel_new)
    return pos_new, vel_new, disp_new, dist_new

#####
coordinates = initial_pos.copy()
velocities = F1.initial_velocities(initial_pos,N,T)
#tables required to compute quantities like forces, energies
displacements = F1.displacement_table(initial_pos[:,[0,1,2]], L)
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
            velocities,timestep, displacements, distances, cutoff,\
            L,K_vecs,Efield)
    
    PE.append(F1.potentialLJ(coordinates,cutoff,L)+F1.Ewald_pot(coordinates,F1.my_legal_kvecs(N_max,L),alpha,V,L))
    PE_EW.append(F1.Ewald_pot(coordinates,F1.my_legal_kvecs(N_max,L),alpha,V,L))
    KE.append(F1.kinetic(coordinates,velocities)[0])
    E.append(F1.potentialLJ(coordinates,cutoff,L) + F1.kinetic(coordinates,velocities)[0])
    P.append(F1.kinetic(coordinates,velocities)[1])
    Temp.append(F1.kinetic(coordinates,velocities)[2])
    
    
    S_k.append(F1.my_calc_sk(F1.my_legal_kvecs(N_max,L),coordinates))
    
    veloci.append(velocities)
    
    G_r.append(F1.my_pair_correlation(F1.dist_ravel(distances),N,num_bins,dr_,L)[0])
    r_.append(F1.my_pair_correlation(F1.dist_ravel(distances),N,num_bins,dr_,L)[1])
    
t1 = time.time()
print(t1-t0)

vacf_arr = np.array([])
num_t = 0
for i in range(len(veloci)):
    add_vel = F1.my_calc_vacf0(veloci,i)
    vacf_arr = np.append(vacf_arr,add_vel)
    num_t=num_t+1
    

#file_KE_write = open('KE_30','w')
#for i in KE:
#    file_KE_write.write(str(i)+'\n')
#    #file_KE_write.write('/')
#file_KE_write.close()
#
#file_PE_write = open('PE_30','w')
#for i in PE:
#    file_PE_write.write(str(i)+'\n')
#    #file_KE_write.write('/')
#file_PE_write.close()
#
#file_PEEW_write = open('PEEW_30','w')
#for i in PE_EW:
#    file_PEEW_write.write(str(i)+'\n')
#    #file_KE_write.write('/')
#file_PEEW_write.close()
#
#file_T_write = open('T_30','w')
#for i in Temp:
#    file_T_write.write(str(i)+'\n')
#    #file_KE_write.write('/')
#file_T_write.close()
#
#file_velo_write = open('velo_30','w')
#for i in veloci:
#    file_velo_write.write(str(i)+'\n')
#    #file_KE_write.write('/')
#file_velo_write.close()
file_velo_write = open('velo_30','w')
for i in vacf_arr:
    file_velo_write.write(str(i)+'\n')
    #file_KE_write.write('/')
file_velo_write.close()