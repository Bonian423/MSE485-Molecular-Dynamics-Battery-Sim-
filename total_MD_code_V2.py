# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 17:21:52 2021

@author: chari
"""

#PYthon functions to import 
#import libraries
import numpy as np
import math
import Lattice_CoO2 as coor

posi = coor.a
T = 300
N = 39

#
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
        qi = i[4]*(1.6e-19)
        for j in pos:
            rj_x = j[0]
            rj_y = j[1]
            rj_z = j[2]
            qj = j[4]*(1.6e-19)
        
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
        qi = i[3]*(1.6e-19)
        V_add_int = (alpha_factor)*(qi)
        V_total_int = V_total_int + V_add_int
    return V_total_int
        
        
def Coulmb_realspace(pos,alph,lbox):
    V_total_real = 0.0
    for i in pos:
        ri_x = i[0]
        ri_y = i[1]
        ri_z = i[2]
        qi = i[4]*(1.6e-19)
        for j in pos:
            rj_x = j[0]
            rj_y = j[1]
            rj_z = j[2]
            qj = j[4]*(1.6e-19)
        
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
    Total_pot = Total_pot*(4.359e-18)
    return Total_pot

def Ewald_force_k(pos,kvecs,alph,vol,lbox):
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-4])
    num_i = 0
    for i in pos:
        qi = i[4]*(1.6e-19)
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = minimum_image(i[0] - j[0],lbox)
            rij_y = minimum_image(i[1] - j[1],lbox)
            rij_z = minimum_image(i[2] - j[2],lbox)
            qj = j[4]*(1.6e-19)
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
    Ewald_force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-4])
    num_i = 0
    for i in pos:
        qi = i[4]*(1.6e-19)
        FX_add = 0.0
        FY_add = 0.0
        FZ_add = 0.0
        for j in pos:
            rij_x = minimum_image(i[0] - j[0],lbox)
            rij_y = minimum_image(i[1] - j[1],lbox)
            rij_z = minimum_image(i[2] - j[2],lbox)
            qj = j[4]*(1.6e-19)
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
    tot_force_arr = tot_force_arr*(8.2387e-3)
    return tot_force_arr

#Functions from HW3/4

def get_temperature(pos,velocities):
    Kb = 1.38e-23
    velo_mag = np.array([])
    
    num = 0
    for i in pos:
        m = i[3]*(1/1000)
        vx = velocities[num,0]
        vy = velocities[num,1]
        vz = velocities[num,2]
        V_mag_sq = vx**2 + vy**2 + vz**2
        mv = m*V_mag_sq
        velo_mag = np.append(velo_mag,mv)
    mean_MV = np.mean(velo_mag)
    T = (mean_MV)/Kb
    return T

def initial_velocities(pos,N,T):
    init_v = np.zeros([N,3])
    Kb = 1.38e-23
    
    v_num = 0
    for i in pos:
        m = i[3]*(1/1000)
        for j in range(0,3):
            random_T1 = np.random.normal()
            init_vel = (np.sqrt((Kb*T)/(m)))*(random_T1)
            init_v[v_num,j] = 0.0
        v_num = v_num + 1
    return init_v
#displacement tables
# no need to change this one from before
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
def kinetic(pos,vel):
    N = np.shape(pos)[0]
    KE_total = 0.0
    P_total = 0.0
    num = 0
    for i in pos:
        m = i[3]*(1/1000)
        vx = vel[num,0]
        vy = vel[num,1]
        vz = vel[num,2]
        V2 = vx**2 + vy**2 + vz**2
        KE_add = (1/2)*(m)*(V2)
        P_add = np.sqrt(2*m*KE_add)
        P_total = P_total + P_add
        KE_total = KE_total + KE_add 
    Temp = (2*KE_total)/(3*N)
    return KE_total,P_total,Temp

## new force and potential for LJ
def potentialLJ(pos,rc,L):
    V_total = 0.0
    for i in pos:
        rx_i = i[0]
        ry_i = i[1]
        rz_i = i[2]
        sigma_i = i[6]
        eps_i = i[5]
        for j in pos:
            rx_j = j[0]
            ry_j = j[1]
            rz_j = j[2]
            sigma_j = j[6]
            eps_j = j[5]
            
            rij_x = minimum_image(rx_i - rx_j,L)
            rij_y = minimum_image(ry_i - ry_j,L)
            rij_z = minimum_image(rz_i - rz_j,L)
            
            eps_ij = (eps_i+eps_j)/2
            sigma_ij = (sigma_i+sigma_j)/2
            
            rij = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
            #print(rij)
            
            if rij >= rc:
                V_add = 0.0
            elif rij == 0.0:
                V_add = 0.0
            else:
                #print(rij)
                V_rc = (4*eps_ij)*((sigma_ij/rc)**6)*((sigma_ij/rc)**6 - 1)
                V_r = (4*eps_ij)*((sigma_ij/rij)**6)*((sigma_ij/rij)**6 - 1) 
                V_add = V_r - V_rc
            
            #print(V_add)
            V_total = V_total + V_add
            #print(V_total)
                
    V_LJ = V_total/2
    V_LJ = V_LJ*(4.359e-18)
    return V_LJ

def forceLJ(pos,rc,L):
    Force_arr = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-4])
    fnum = 0
    for i in pos:
        rx_i = i[0]
        ry_i = i[1]
        rz_i = i[2]
        sigma_i = i[6]
        eps_i = i[5]
        
        Fx = 0.0
        Fy = 0.0
        Fz = 0.0
        for j in pos:
            rx_j = j[0]
            ry_j = j[1]
            rz_j = j[2]
            sigma_j = j[6]
            eps_j = j[5]
            
            rij_x = minimum_image(rx_i-rx_j,L)
            rij_y = minimum_image(ry_i-ry_j,L)
            rij_z = minimum_image(rz_i-rz_j,L)
            
            eps_ij = (eps_i+eps_j)/2
            sigma_ij = (sigma_i+sigma_j)/2
            
            rij = np.sqrt(rij_x**2 + rij_y**2 + rij_z**2)
            if rij >= rc or rij == 0.0:
                F_mag = 0.0
            else:
                F_mag = (24*eps_ij)*(sigma_ij/(rij**2))*((sigma_ij/rij)**6)*(2*((sigma_ij/rij)**6)-1)
            
            Fx_add = (F_mag)*(rij_x)
            Fy_add = (F_mag)*(rij_y)
            Fz_add = (F_mag)*(rij_z)
            
            Fx = Fx + Fx_add
            Fy = Fy + Fy_add
            Fz = Fz + Fz_add
            
        Force_arr[fnum,0] = Fx
        Force_arr[fnum,1] = Fy
        Force_arr[fnum,2] = Fz
        fnum = fnum + 1
    Force_arr = Force_arr*(8.2387e-8)
    return Force_arr

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

def external_force(pos,E):
    F_ext = np.zeros([np.shape(pos)[0],np.shape(pos)[1]-4])
    fnum = 0
    for i in pos:
        q = i[4]*(1.6e-19)
        Fq = q*E
        F_ext[fnum] = Fq
        fnum = fnum + 1
    return F_ext