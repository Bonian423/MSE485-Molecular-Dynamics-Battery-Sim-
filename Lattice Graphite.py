# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 19:40:48 2021

@author: Bonia
"""
# -*- coding: utf-8 -*-

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

"""lattice parameters has format [x,y,z]"""
"""Unit of the lattice parameters: Angstrom"""
lattice_params = np.array([2.46772,2.46772*np.sin(np.pi/3), 8.68504])

"""structure parameters has format [x,y,z]"""
"""structure parameters are formatted as ratios of lattice parameter"""
C_struct_param = np.array([[0.00,    0.00,    0.25 ],[0.0,    0.0,    0.75],[0.66667,    0.33333,    0.25 ],[0.33333,    0.66667,    0.75]])

"""Atom Positions in unit cell in angstrom"""
C_UC_pos = np.multiply(C_struct_param, lattice_params)

for atom in C_UC_pos:
    x0 = atom[0]
    y0 = atom[1]
    """this conversion is proved to be unnecessary"""
    atom[0] = x0-y0*np.cos(np.pi/3)
    atom[1] = y0*np.sin(np.pi/3)


"""properties of atoms at each position"""

Li_col = np.array([0.1,0.1])
Co_col = np.array([0.9,0.9])
O_col = np.array([0.5,0.5,0.5,0.5])
fixed_UC_color = np.concatenate((Co_col, O_col),axis = 0)

def lattice_construction(xnum, ynum, znum):
    """
    Construct a fixed lattice for SiO2 by defining how many unit cells we want in each dimension

    Parameters
    ----------
    xnum : int, >=1
        number of unit cells in the x direction
    ynum : int, >=1
        number of unit cells in the y direction
    znum : int, >=1
        number of unit cells in the z direction
    Returns
    -------
    np.array structure with all the atoms in the constructed lattice
    """
    xlat = C_UC_pos
    for x in range(1,xnum):
        xshift = lattice_params[0]
        shift = np.array([x*xshift,0,0])
        xlat = np.concatenate((xlat,np.add(shift, C_UC_pos)),axis =0)
    ylat = xlat    
    print(xlat.shape)
    for y in range(1,ynum):
        yshift = lattice_params[1]*np.sin(np.pi/3)
        xshift = 2.93855*np.cos(np.pi/3)
        shift = np.array([y*xshift,y*yshift,0])
        ylat = np.concatenate((ylat,np.add(shift, xlat)),axis =0)
    print(ylat.shape)
    zlat = ylat
    for z in range(1,znum):
        zshift = lattice_params[2]
        shift = np.array([0,0,z*zshift])
        zlat = np.concatenate((zlat,np.add(shift, ylat)),axis =0)        
    print(zlat.shape)
    return zlat

def minimum_image(r, L):
  return r - L*np.round(r / L)



dim = 15
a = lattice_construction(7, 8, 2)
a = minimum_image(a, dim)


plt.figure(1)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2])
ax.set_xlim3d(left=-dim/2, right=dim/2)
ax.set_ylim3d(bottom=-dim/2, top=dim/2)
ax.set_zlim3d(bottom=-dim/2, top=dim/2)
ax.view_init(-180, 90)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2])
ax.set_xlim3d(left=-dim/2, right=dim/2)
ax.set_ylim3d(bottom=-dim/2, top=dim/2)
ax.set_zlim3d(bottom=-dim/2, top=dim/2)
ax.view_init(90, 270)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.figure(3)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2])
ax.set_xlim3d(left=-dim/2, right=dim/2)
ax.set_ylim3d(bottom=-dim/2, top=dim/2)
ax.set_zlim3d(bottom=-dim/2, top=dim/2)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

"""Unit cell visualization"""
# plt.figure(3)
# bx = plt.axes(projection='3d')
# bx.scatter3D(fixed_UC_pos[:,0],fixed_UC_pos[:,1],fixed_UC_pos[:,2], c = fixed_UC_color)
# bx.set_xlim3d(left=0, right=10)
# bx.set_ylim3d(bottom=0, top=10)
# bx.set_zlim3d(bottom=0, top=10)
# bx.set_xlabel("x")
# bx.set_ylabel("y")
# bx.set_zlabel("z")
# bx.view_init(-90, 90)
# plt.figure(4)
# bx = plt.axes(projection='3d')
# bx.scatter3D(fixed_UC_pos[:,0],fixed_UC_pos[:,1],fixed_UC_pos[:,2], c = fixed_UC_color)
# bx.set_xlim3d(left=0, right=10)
# bx.set_ylim3d(bottom=0, top=10)
# bx.set_zlim3d(bottom=0, top=10)
# bx.set_xlabel("x")
# bx.set_ylabel("y")
# bx.set_zlabel("z")
