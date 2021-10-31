# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 15:26:12 2021

@author: Bonia
"""
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

"""lattice parameters has format [x,y,z]"""
"""Unit of the lattice parameters: Angstrom"""
lattice_params = np.array([2.93855,2.93855, 9.85993])

"""structure parameters has format [x,y,z]"""
"""structure parameters are formatted as ratios of lattice parameter"""
Li_struct_param = np.array([[0.66667,    0.33333,    0.25821 ],[0.33333,    0.66667,    0.75821]])
Co_struct_param = np.array([[0.66667,    0.33333,    0.00051 ],[0.33333,    0.66667,    0.50051]])
O_struct_param = np.array([[0.00000,    0.00000,    0.88677],[0.00000,    0.00000,    0.38677],[0.66667,    0.33333,    0.61344],[0.33333,    0.66667,    0.11344]])

"""Atom Positions in unit cell in angstrom"""
Li_UC_pos = np.multiply(Li_struct_param, lattice_params)
Co_UC_pos = np.multiply(Co_struct_param, lattice_params)
O_UC_pos = np.multiply(O_struct_param, lattice_params)

fixed_UC_pos = np.concatenate((Co_UC_pos, O_UC_pos),axis = 0)
for atom in fixed_UC_pos:
    x0 = atom[0]
    y0 = atom[1]
    atom[0] = x0+y0/np.tan(np.pi/3)
    atom[1] = y0*np.sin(np.pi/3)


"""properties of atoms at each position"""
Li_properties_list = np.array([["AASHKAJ","sjadahkhd"],["sjshadjkaj","ajksdhak"]])
Co_properties_list = np.array([["AASHKAJ","sjadahkhd"],["sjshadjkaj","ajksdhak"]])
O_properties_list = np.array([["AASHKAJ","sjadahkhd"],["sjshadjkaj","ajksdhak"],["sjshadjkaj","ajksdhak"],["sjshadjkaj","ajksdhak"]])

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
    xlat = fixed_UC_pos
    for x in range(1,xnum):
        xshift = lattice_params[0]
        shift = np.array([x*xshift,0,0])
        xlat = np.concatenate((xlat,np.add(shift, fixed_UC_pos)),axis =0)
    ylat = xlat    
    print(xlat.shape)
    for y in range(1,ynum):
        yshift = lattice_params[1]
        shift = np.array([0,y*yshift,0])
        ylat = np.concatenate((ylat,np.add(shift, xlat)),axis =0)
    print(ylat.shape)
    zlat = ylat    
    for z in range(1,znum):
        zshift = lattice_params[2]
        shift = np.array([0,0,z*zshift])
        zlat = np.concatenate((zlat,np.add(shift, ylat)),axis =0)
    print(zlat.shape)
    return zlat

a = lattice_construction(5, 5, 2)
plt.figure(1)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2])
plt.figure(2)
bx = plt.axes(projection='3d')
bx.scatter3D(fixed_UC_pos[:,0],fixed_UC_pos[:,1],fixed_UC_pos[:,2])
bx.set_xlim3d(left=0, right=10)
bx.set_ylim3d(left=0, right=10)
bx.set_zlim3d(left=0, right=10)
bx.set_xlabel("x")
bx.set_ylabel("y")
bx.set_zlabel("z")