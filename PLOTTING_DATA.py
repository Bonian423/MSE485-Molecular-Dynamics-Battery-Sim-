# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:00:43 2021

@author: chari
"""

#Plotting 
import numpy as np
import matplotlib.pyplot as plt
import total_MD_code_V2 as F1

time_steps = np.arange(0,100)

KE_arr = np.array([])
PE_arr = np.array([])
PEEW_arr = np.array([])
T_arr = np.array([])
vacf_0 = np.array([])

with open('KE_30') as f:
    for line in f:
        KE_arr = np.append(KE_arr,float(line))
        
with open('PE_30') as f:
    for line in f:
        PE_arr = np.append(PE_arr,float(line))

with open('PEEW_30') as f:
    for line in f:
        PEEW_arr = np.append(PEEW_arr,float(line))

with open('T_30') as f:
    for line in f:
        T_arr = np.append(T_arr,float(line))        

with open('velo_30') as f:
    for line in f:
        vacf_0 = np.append(vacf_0,float(line))

plt.figure(1)
plt.plot(np.arange(0,len(KE_arr)),KE_arr,'--bo')
plt.title('Kinetic Energy vs Time Step')
plt.ylabel('KE (J)')
plt.xlabel('Time Step (1e-12)')
        
        
plt.figure(2)
plt.plot(np.arange(0,len(PE_arr)),PE_arr,'--bo')
plt.title('Total Potential Energy vs Time Step')
plt.ylabel('PE (J)')
plt.xlabel('Time Step (1e-12)')
                
plt.figure(3)
plt.plot(np.arange(0,len(PEEW_arr)),PEEW_arr,'--bo')
plt.title('Ewald Potential Energy vs Time Step')
plt.ylabel('PE (J)')
plt.xlabel('Time Step (1e-12)')        
        
plt.figure(4)
plt.plot(np.arange(0,len(T_arr)),T_arr,'--bo')
plt.title('MD Temperature vs Time Step')
plt.ylabel('Temperature (K)')
plt.xlabel('Time Step (1e-12)') 

plt.figure(6)
plt.plot(np.arange(0,len(vacf_0)),vacf_0,'--bo',label = 'D='+str(np.round(F1.my_diffusion_constant(vacf_0),3)))
plt.title('V-V autocorrelation function')
plt.xlabel('Time Step (1e-12)') 
plt.legend()
            
       
        
        
#Gr_avg = np.array([])
#Sk_avg = np.array([])
#for i in S_k:
#    add_avgS = np.average(i)
#    Sk_avg = np.append(Sk_avg,add_avgS)
#for i in G_r:
#    add_avgG = np.average(i)
#    Gr_avg = np.append(Gr_avg,add_avgG)
     
        
        
        
###
import matplotlib.pyplot as plt
import numpy as np
for i in E:
    if i > 0.1:
        E.remove(i)
plt.plot(np.arange(0,len(E)),E)  
plt.title('Total Energy vs Time Step - T1')
plt.ylabel('Energy (J)')
plt.xlabel('Time Step (1e-12)')      
plt.savefig('spydata/Trial1_E',dpi=500,bbox_inches = "tight") 
#        
        
        
        
        