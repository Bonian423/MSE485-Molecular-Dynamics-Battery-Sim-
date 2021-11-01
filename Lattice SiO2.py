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
    """this conversion is proved to be unnecessary"""
    atom[0] = x0-y0*np.cos(np.pi/3)
    atom[1] = y0*np.sin(np.pi/3)


"""properties of atoms at each position"""
Li_properties_list = np.array([["G"],["G"]])
Co_properties_list = np.array([["R"],["R"]])
O_properties_list = np.array([["B"],["B"],["B"],["B"]])
fixed_UC_prop = np.concatenate((Co_UC_pos, O_UC_pos),axis = 0)
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
    xlat = fixed_UC_pos
    xcol = fixed_UC_color
    for x in range(1,xnum):
        xshift = lattice_params[0]
        shift = np.array([x*xshift,0,0])
        xlat = np.concatenate((xlat,np.add(shift, fixed_UC_pos)),axis =0)
        xcol = np.concatenate((xcol, fixed_UC_color),axis =0)
    ylat = xlat    
    ycol = xcol
    print(xlat.shape)
    for y in range(1,ynum):
        yshift = lattice_params[1]*np.sin(np.pi/3)
        xshift = 2.93855*np.cos(np.pi/3)
        shift = np.array([y*xshift,y*yshift,0])
        ylat = np.concatenate((ylat,np.add(shift, xlat)),axis =0)
        ycol = np.concatenate((xcol, ycol),axis =0)     
    print(ylat.shape)
    zlat = ylat
    zcol = ycol    
    for z in range(1,znum):
        zshift = lattice_params[2]
        shift = np.array([0,0,z*zshift])
        zlat = np.concatenate((zlat,np.add(shift, ylat)),axis =0)        
        zcol = np.concatenate((zcol, ycol),axis =0)        
    print(zlat.shape)
    return zlat, zcol
def minimum_image(r, L):
    PBC_pos = np.zeros([len(r),3])

    num_i = 0
    for i in r:
        num_j = 0
        for j in i:
            if -L/2 <= j < L/2:
                PBC_pos[num_i,num_j] = j
            if j >= L/2:
                j_add1 = j - (L)*(j//(L)+1)
                if -L/2 <= j_add1 < L/2:
                    PBC_pos[num_i,num_j] = j_add1
                else:
                    PBC_pos[num_i,num_j]  = j_add1 + L
            if j < -L/2:
                j_add1 = j - (L)*(j//(L)+1)
                if -L/2 <= j_add1 < L/2:
                    PBC_pos[num_i,num_j] = j_add1
                else:
                    PBC_pos[num_i,num_j]  = j_add1 + L    
    
            num_j = num_j + 1
        num_i = num_i+1
    
    return PBC_pos
    pass

dim = 15
a , b = lattice_construction(5, 6, 2)
a = minimum_image(a, dim)
plt.figure(1)
dim = 15
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2], c = b)
ax.set_xlim3d(left=-dim/2, right=dim/2)
ax.set_ylim3d(bottom=-dim/2, top=dim/2)
ax.set_zlim3d(bottom=-dim/2, top=dim/2)
ax.view_init(-180, 90)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2], c = b)
ax.set_xlim3d(left=-dim/2, right=dim/2)
ax.set_ylim3d(bottom=-dim/2, top=dim/2)
ax.set_zlim3d(bottom=-dim/2, top=dim/2)
ax.view_init(90, 270)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.figure(3)
ax = plt.axes(projection='3d')
ax.scatter3D(a[:,0], a[:,1], a[:,2], c = b)
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
