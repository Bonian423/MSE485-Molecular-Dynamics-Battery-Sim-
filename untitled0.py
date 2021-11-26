# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 19:08:11 2021

@author: Bonia
"""
    for x in range(0,pos_new.shape[0]):
        if pos_new[x][4] == 1:  
            pos_new[x][0:3]= pos[x,[0,1,2]] + dt*vel_half[x]
            
            
            
