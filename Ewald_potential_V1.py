# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 22:22:05 2021

@author: chari
"""

#Ewald test
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import cmath
#
def cubic_lattice(tiling, L):
    x = np.linspace(-0.35*L,0.35*L,tiling)
    y = np.linspace(-0.35*L,0.35*L,tiling)
    z = np.linspace(-0.35*L,0.35*L,tiling)
    
    x_n,y_n,z_n = np.meshgrid(x,y,z,indexing='ij')
    
    lat = np.column_stack((x_n.ravel(),y_n.ravel(),z_n.ravel()))
    
    return lat
    pass
def pos_in_box(pos, lbox):
  return (pos+lbox/2.) % lbox - lbox/2.
def minimum_image(r, L):
  return r - L*np.round(r / L)
def my_legal_kvecs(maxn, lbox):
    k_fact = (2*np.pi)/(lbox)

    x = np.linspace(0.,maxn,maxn+1)
    kx_ = (k_fact)*(x)

    y = np.linspace(0.,maxn,maxn+1)
    ky_ = (k_fact)*(y)

    z = np.linspace(0.,maxn,maxn+1)
    kz_ = (k_fact)*(z)

    kx, ky, kz = np.meshgrid(kx_, ky_, kz_, indexing='ij')
    
    kvecs = np.column_stack((kx.ravel(),ky.ravel(),kz.ravel()))
    
    return kvecs

###############################################################################



#
lbox = 5
#alpha = 100/(lbox**2)
cutoff = lbox/2


N = 8
V = lbox**3
alpha = (np.pi)*(N/(V**2))**(1/3)


#
init_pos = cubic_lattice(2,lbox)
Qs = np.array([[1],[-1],[1],[1],[-1],[1],[-1],[1]])
init_pos = np.append(init_pos,Qs,axis=1)
init_kvecs = my_legal_kvecs(4,lbox)

fig1 = plt.figure(2)
ax1 = plt.axes(projection='3d')
scatP1 = ax1.scatter3D(init_pos[:,0],init_pos[:,1],init_pos[:,2])
plt.show()
#

def displacement_table(coordinates, L):
    coor_shape1 = np.shape(coordinates)[0]
    
    disp_table = np.zeros([coor_shape1,coor_shape1,3])
    num1 = 0
    for i in coordinates:
        x_i = i[0]
        y_i = i[1]
        z_i = i[2]
        num2 = 0
        for j in coordinates:
            x_j = j[0]
            y_j = j[1]
            z_j = j[2]
        
            xij = x_i - x_j
            yij = y_i - y_j
            zij = z_i - z_j
            
            pre_MIC = np.array([xij,yij,zij])
            
            post_MIC = pos_in_box(pre_MIC,L)
            
            disp_table[num1,num2,0] = post_MIC[0]
            disp_table[num1,num2,1] = post_MIC[1]
            disp_table[num1,num2,2] = post_MIC[2]
            
            num2 = num2 + 1
        num1 = num1 + 1
        
    return disp_table
    pass

init_disp = displacement_table(init_pos,lbox)


#            
            
        
def Coulmb_kspace(pos,kvecs,alph,vol):
    V_total = 0.0
    for i in pos:
        ri_x = i[0]
        ri_y = i[1]
        ri_z = i[2]
        qi = i[3]
        for j in pos:
            rj_x = j[0]
            rj_y = j[1]
            rj_z = j[2]
            qj = j[3]
        
            rij_x = minimum_image(ri_x - rj_x,lbox)
            rij_y = minimum_image(ri_y - rj_y,lbox)
            rij_z = minimum_image(ri_z - rj_z,lbox)
        
            rij_mag = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
        
            if rij_mag == 0.0:
                continue
            else:
                for k in kvecs:
                    kx = k[0]
                    ky = k[1]
                    kz = k[2]
                
                    dot_KR = (kx*rij_x) + (ky*rij_y) + (kz*rij_z)
                    k_mag = np.sqrt(kx**2 + ky**2 + kz**2)
                    #print(k_mag)
                    if k_mag == 0.0:
                        continue
                    else:
                        V_add = (1/(2*(vol)))*(qi*qj)*((4*np.pi)/(k_mag**2))*(np.exp(-(k_mag**2)/(4*alph)))*(abs(complex(math.cos(dot_KR),math.sin(dot_KR))))**2
                        #print(V_add)
                        V_total = V_total + V_add
    return V_total

def Coulmb_selfInteract(pos,alph):
    V_total_int = 0.0
    alpha_factor = (-alph)/(np.sqrt(np.pi))
    for i in pos:
        qi = i[3]
        V_add_int = (alpha_factor)*(qi)
        V_total_int = V_total_int + V_add_int
    return V_total_int
        
        
def Coulmb_realspace(pos,alph):
    V_total_real = 0.0
    for i in pos:
        ri_x = i[0]
        ri_y = i[1]
        ri_z = i[2]
        qi = i[3]
        for j in pos:
            rj_x = j[0]
            rj_y = j[1]
            rj_z = j[2]
            qj = j[3]
        
            rij_x = minimum_image(ri_x - rj_x,lbox)
            rij_y = minimum_image(ri_y - rj_y,lbox)
            rij_z = minimum_image(ri_z - rj_z,lbox)
        
            rij_mag = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
            if rij_mag == 0.0:
                continue
            else:
                V_add_real = (1/2)*((qi*qj)/(rij_mag))*math.erfc(np.sqrt(alph)*(rij_mag))
                V_total_real = V_total_real + V_add_real
    return V_total_real
    

def Ewald_pot(pos,kvecs,alph,vol):
    Total_pot = Coulmb_kspace(pos,kvecs,alph,vol) - Coulmb_selfInteract(pos,alph) + Coulmb_realspace(pos,alph)
    return Total_pot


# Forces
    
def Ewald_force_k(pos,kvecs,alph,vol):
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-1])
    num_i = 0
    for i in pos:
        qi = i[3]
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = i[0] - j[0]
            rij_y = i[1] - j[1]
            rij_z = i[2] - j[2]
            qj = j[3]
            rij_mag = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
            if rij_mag == 0.0:
                continue
            else:
                for k in kvecs:
                    kx = k[0]
                    ky = k[1]
                    kz = k[2]
                    dot_KR = (kx*rij_x) + (ky*rij_y) + (kz*rij_z)
                    k_mag = np.sqrt(kx**2 + ky**2 + kz**2)
                    if k_mag == 0.0:
                        continue
                    else:
                        fx_add = (qi*qj)*(1/vol)*(4*np.pi*kx)*(1/(k_mag**2))*(np.exp((-k_mag**2)/(4*alph)))*(np.sin(dot_KR))
                        fy_add = (qi*qj)*(1/vol)*(4*np.pi*ky)*(1/(k_mag**2))*(np.exp((-k_mag**2)/(4*alph)))*(np.sin(dot_KR))
                        fz_add = (qi*qj)*(1/vol)*(4*np.pi*kz)*(1/(k_mag**2))*(np.exp((-k_mag**2)/(4*alph)))*(np.sin(dot_KR))
                        
                        FX_add = FX_add + fx_add
                        FY_add = FY_add + fy_add
                        FZ_add = FZ_add + fz_add
        Ewald_force_arr[num_i][0] = FX_add
        Ewald_force_arr[num_i][1] = FY_add
        Ewald_force_arr[num_i][2] = FY_add
        num_i = num_i+1
    
    return Ewald_force_arr


def Ewald_force_r(pos,alph):
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-1])
    num_i = 0
    for i in pos:
        qi = i[3]
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = i[0] - j[0]
            rij_y = i[1] - j[1]
            rij_z = i[2] - j[2]
            qj = j[3]
            rij_mag = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
            if rij_mag == 0.0:
                continue
            else:
                fx_add = (qi*qj)*((2*np.sqrt((alph)/(np.pi)))*(np.exp(-alph*(rij_mag**2))) + (1/rij_mag)*(math.erfc(np.sqrt(alph)*(rij_mag))))*(rij_x/(rij_mag**2))
                fy_add = (qi*qj)*((2*np.sqrt((alph)/(np.pi)))*(np.exp(-alph*(rij_mag**2))) + (1/rij_mag)*(math.erfc(np.sqrt(alph)*(rij_mag))))*(rij_y/(rij_mag**2))
                fz_add = (qi*qj)*((2*np.sqrt((alph)/(np.pi)))*(np.exp(-alph*(rij_mag**2))) + (1/rij_mag)*(math.erfc(np.sqrt(alph)*(rij_mag))))*(rij_z/(rij_mag**2))
                
                FX_add = FX_add + fx_add
                FY_add = FY_add + fy_add
                FZ_add = FZ_add + fz_add
        Ewald_force_arr[num_i][0] = FX_add
        Ewald_force_arr[num_i][1] = FY_add
        Ewald_force_arr[num_i][2] = FZ_add
        num_i = num_i + 1
    return Ewald_force_arr
    
def Ewald_force(pos,kvecs,alph,vol):
    tot_force_arr = Ewald_force_r(pos,alph) + Ewald_force_k(pos,kvecs,alph,vol)
    return tot_force_arr

