# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 14:50:59 2021

@author: chari
"""

#Project MD total function
#import libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
import cmath

#import pythonn functions needed
#Ewald necessary functions
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

def Coulmb_kspace(pos,kvecs,alph,vol,lbox):
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
        
        
def Coulmb_realspace(pos,alph,lbox):
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
    

def Ewald_pot(pos,kvecs,alph,vol,lbox):
    Total_pot = Coulmb_kspace(pos,kvecs,alph,vol,lbox) - Coulmb_selfInteract(pos,alph) + Coulmb_realspace(pos,alph,lbox)
    return Total_pot

def Ewald_force_k(pos,kvecs,alph,vol,lbox):
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-1])
    num_i = 0
    for i in pos:
        qi = i[3]
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = minimum_image(i[0] - j[0],lbox)
            rij_y = minimum_image(i[1] - j[1],lbox)
            rij_z = minimum_image(i[2] - j[2],lbox)
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


def Ewald_force_r(pos,alph,lbox):
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-1])
    num_i = 0
    for i in pos:
        qi = i[3]
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = minimum_image(i[0] - j[0],lbox)
            rij_y = minimum_image(i[1] - j[1],lbox)
            rij_z = minimum_image(i[2] - j[2],lbox)
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
    
def Ewald_force(pos,kvecs,alph,vol,lbox):
    tot_force_arr = Ewald_force_r(pos,alph,lbox) + Ewald_force_k(pos,kvecs,alph,vol,lbox)
    return tot_force_arr

#Functions from HW3/4
#Temperature/init velocity
def get_temperature(mass, velocities):
    Kb = 1.0
    velo_mag = np.array([])
    for i in velocities:
        vx = i[0]
        vy = i[1]
        vz = i[2]
    
        V_mag = vx**2 + vy**2 + vz**2
    
        velo_mag_2 = np.append(velo_mag,V_mag)
    
    Vsq_mean = np.mean(velo_mag_2)
    T = ((mass)*(Vsq_mean))/Kb
    
    return T
    pass

def initial_velocities(N, m, T):  
    init_v = np.zeros([N,3])
    Kb = 1.0
    
    vnum1 = 0
    for i in init_v:
        vnum2 = 0
        for j in i:
            random_T1 = np.random.normal()
            MB_init_v = (np.sqrt((Kb*T)/(m)))*(random_T1)
            init_v[vnum1,vnum2] = MB_init_v
            vnum2 = vnum2+1
        vnum1 = vnum1+1
    
    return init_v
    pass

#displacement tables
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

#Kinetic Energy
def kinetic(m, v, N):
    KE_array =  np.zeros([np.shape(v)[0]])
    
    KE_count = 0
    for i in v:
        vx = i[0]
        vy = i[1]
        vz = i[2]
        Vel_tot = np.sqrt(vx**2 + vy**2 + vz**2)
        KE_add = (1/2)*(m)*((Vel_tot)**2)
        
        KE_array[KE_count] = KE_add
        
        KE_count = KE_count + 1
    
    KE_total = np.sum(KE_array)
    P_total = np.sqrt(2*m*KE_total)
    Temp = (2*KE_total)/(3*N)
    
    return KE_total,P_total,Temp
    pass

#Lennard Jones Force
def potential(dist, rc):
    eps = 1.0
    sigma = 1.0
    V_array = np.array([])
    for i in dist:
        for j in i:
            if j == 0.0:
                continue
            elif j >= rc:
                V_add = 0.0
                V_array = np.append(V_array,V_add)
            elif j < rc:
                V_rc = (4*eps)*((sigma/rc)**6)*((sigma/rc)**6 - 1)
                V_r = (4*eps)*((sigma/j)**6)*((sigma/j)**6 - 1) 
                V_add = V_r - V_rc
                V_array = np.append(V_array,V_add)

    V_sum_all = np.sum(V_array)

    V_total = V_sum_all/2
    
    return V_total
    pass

def force(disp, dist, rc):
    eps = 1.0
    sigma = 1.0
    
    Force_array = np.zeros([np.shape(disp)[0],3])

    fnum_1 = 0
    for i in disp:
        F_i_x_add = 0.0
        F_i_y_add = 0.0
        F_i_z_add = 0.0
        for j in i:
            rx = j[0]
            ry = j[1]
            rz = j[2]
    
            r_mag = np.sqrt(rx**2 + ry**2 + rz**2)
        
            if r_mag > rc or r_mag == 0.0:
                F_mag = 0.0
            else:
                F_mag = (24*eps)*(sigma/(r_mag**2))*((sigma/r_mag)**6)*(2*((sigma/r_mag)**6)-1)
        
            F_j_x_add = (F_mag)*(rx)
            F_j_y_add = (F_mag)*(ry)
            F_j_z_add = (F_mag)*(rz)
        
            F_i_x_add = F_i_x_add + F_j_x_add
            F_i_y_add = F_i_y_add + F_j_y_add
            F_i_z_add = F_i_z_add + F_j_z_add
        
        Force_array[fnum_1,0] = F_i_x_add
        Force_array[fnum_1,1] = F_i_y_add
        Force_array[fnum_1,2] = F_i_z_add
    
        fnum_1 = fnum_1 + 1
    return Force_array
    pass
## Extra calculations
def my_pair_correlation(dists, natom, nbins, dr, lbox):
    histogram = np.histogram(dists, bins=nbins, range=(0, nbins*dr))
    r = (histogram[1] + dr/2)[:-1] # centers of the bins
    
    init_V = lbox**3
    density = natom/init_V
    
    Vol = np.array([])
    for i in r:
        add = ((4*np.pi)/3)*((i + dr/2)**3 - (i - dr/2)**3)
        Vol = np.append(Vol,add)
    
    ideal_gas_arr = np.array([])
    for i in Vol:
        add1 = density*i
        ideal_gas_arr = np.append(ideal_gas_arr,add1)
        
    final = np.array([])
    for i in range(len(Vol)):
        add_g = (histogram[0][i]/ideal_gas_arr[i])/((natom-1)/2)
        final = np.append(final,add_g)
        
    return final,r
    pass

def my_calc_rhok(kvecs, pos):
    rho = np.array([])
    for i in kvecs:
        kx = i[0]
        ky = i[1]
        kz = i[2]
        dot_add = 0
        for j in pos:
            posx = j[0]
            posy = j[1]
            posz = j[2]
        
            add1 = (kx*posx) + (ky*posy) + (kz*posz)
            add = complex(np.cos(add1),(-1)*np.sin(add1))
        
            dot_add = dot_add + add
        rho = np.append(rho,dot_add)
        dot_add = 0
    
    return rho
    pass

def my_calc_sk(kvecs, pos):
    N = len(pos)
    
    rho_k = my_calc_rhok(kvecs,pos)
    rho_min_k =  my_calc_rhok(-kvecs,pos)
    
    S_k = np.array([])
    
    num2 = 0
    for i in rho_k:
        add = (1/N)*np.mean(i*rho_min_k[num2])
        S_k = np.append(S_k,add)
        num2 = num2+1
        
    return S_k
    pass

def my_calc_vacf0(all_vel, t):
    N = len(all_vel[0])
    v_0_vec = all_vel[0]
    v_t_vec = all_vel[t]
    
    Dot_arr = 0
    num_vec = 0
    for i in v_0_vec:
        vx_0 = i[0]
        vy_0 = i[1]
        vz_0 = i[2]
        
        vx_t = v_t_vec[num_vec,0]
        vy_t = v_t_vec[num_vec,1]
        vz_t = v_t_vec[num_vec,2]
        
        add = (vx_0*vx_t) + (vy_0*vy_t) + (vz_0*vz_t)
        
        Dot_arr = Dot_arr + add
        
        num_vec = num_vec+1
    
    vacf0 = (Dot_arr)/N
    
    return vacf0

def my_diffusion_constant(vacf):
    dt = 0.032 #NOTE Do not change this value.
    inte = np.trapz(vacf,dx=dt)
    D = (1/3)*(inte)
    
    return D

def dist_ravel(dist):
    ravel = np.ravel(dist)
    ravel2 = np.array([])
    for i in ravel:
        if i == 0:
            continue
        else:
            ravel2 = np.append(ravel2,i)
    return ravel2



