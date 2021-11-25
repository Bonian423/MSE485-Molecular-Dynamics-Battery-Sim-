# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:16:21 2021

@author: Bonia
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 30 15:26:12 2021

@author: Bonia
"""
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


class lattice: 
    
    
    def lattice_construction(self, numofcells, UC_ary, dim, fixed_UC_color):
        """
        Construct a fixed lattice for SiO2 by defining how many unit cells we want in each dimension
        Parameters
        ----------
        numofcells: list of dimensions
            xnum : int, >=1
                number of unit cells in the x direction
            ynum : int, >=1
                number of unit cells in the y direction
            znum : int, >=1
                number of unit cells in the z direction
        Returns: 
            zlat: final lattice
            zcol: color matrix
            zprop: property matrix of each atom
            * all have the same index for the same atom
        -------
        np.array structure with all the atoms in the constructed lattice
        """

        xnum = numofcells[0]
        ynum = numofcells[1]
        znum = numofcells[2]
        fixed_UC_pos = UC_ary[:,0:3]
        fixed_UC_prop = UC_ary[:,3:]
        xlat = fixed_UC_pos
        xprop = fixed_UC_prop
        xcol = fixed_UC_color
        for x in range(1,xnum):
            xshift = self.lattice_params[0]
            shift = np.array([x*xshift,0,0])
            xlat = np.concatenate((xlat,np.add(shift, fixed_UC_pos)),axis =0)
            xcol = np.concatenate((xcol, fixed_UC_color),axis =0)
            xprop = np.vstack([xprop, fixed_UC_prop])
        ylat = xlat
        yprop = xprop    
        ycol = xcol
        for y in range(1,ynum):
            yshift = self.lattice_params[1]*np.sin(np.pi/3)
            xshift = self.lattice_params[0]*np.cos(np.pi/3)
            shift = np.array([y*xshift,y*yshift,0])
            ylat = np.concatenate((ylat,np.add(shift, xlat)),axis =0)
            ycol = np.concatenate((xcol, ycol),axis =0)     
            yprop = np.vstack([yprop, xprop])
        zlat = ylat
        zcol = ycol    
        zprop = yprop  
        for z in range(1,znum):
            zshift = self.lattice_params[2]
            shift = np.array([0,0,z*zshift])
            zlat = np.concatenate((zlat,np.add(shift, ylat)),axis =0)        
            zcol = np.concatenate((zcol, ycol),axis =0)        
            zprop = np.vstack([zprop, yprop])      
        return np.concatenate((zlat, zprop),axis = 1), zcol
    def shift_to_center(self,r,dim):
        """shift lattice to the center of the space before applying minimum image
        ---------
        input:
            r: lattice
            dim: real space dimension of the entire lattice
        output: shifted lattice      
        """
        xshift = dim[0]/2
        yshift = dim[1]/2
        for rtemp in r:
            rtemp[0] -= xshift
            rtemp[1] -= yshift
        return r
    def minimum_image(self,r, L):
        Lx = L[0]
        Ly = L[1]
        Lz = L[2]
        r[:,0] = r[:,0] - Lx*np.round(r[:,0] / Lx)
        r[:,1] = r[:,1] - Ly*np.round(r[:,1] / Ly)
        r[:,2] = r[:,2] - Lz*np.round(r[:,2] / Lz)
        return r
    
    def crop_cubic(self,r,mindim,col):
        """shift lattice to the center of the space before applying minimum image
        ---------
        input:
            r: lattice
            mindim: minimum dimension of the entire lattice
        output: cropped lattice based on the minimum dimension     
        """
        oob = []
        for i in range(r.shape[0]):
            if np.abs(r[i][0]) > mindim/2 or np.abs(r[i][1]) > mindim/2 or np.abs(r[i][2]) > mindim/2:
                oob.append(i)
        r = np.delete(r,oob,0)
        col = np.delete(col,oob)
        return r,col
    
    def removedup(self,r):
        rnew = r
        dellist = []
        for i in range (0,r.shape[0]):
            for j in range (i+1,r.shape[0]):
                rij = np.abs(np.subtract(r[j],r[i]))
                if np.sum(rij) <  0.0001:
                    dellist.append(j)
        rnew = np.delete(rnew,dellist,0)
        return rnew   

    def init(self, numofcells,boxdim):
        """lattice parameters has format [x,y,z]"""
        """Unit of the lattice parameters: Angstrom"""
        self.lattice_params = np.array([2.46772,2.46772, 8.68504])
    
        """structure parameters has format [x,y,z]"""
        """structure parameters are formatted as ratios of lattice parameter"""
        C_struct_param = np.array([[0.00,    0.00,    0.25 ],[0.0,    0.0,    0.75],[0.66667,    0.33333,    0.25 ],[0.33333,    0.66667,    0.75]])
  
    
        """Atom Positions in unit cell in angstrom"""
        C_UC_pos = np.multiply(C_struct_param, self.lattice_params)

    
        fixed_UC_pos = C_UC_pos
        for atom in fixed_UC_pos:
            x0 = atom[0]
            y0 = atom[1]
            atom[0] = x0-y0*np.cos(np.pi/3)
            atom[1] = y0*np.sin(np.pi/3)
    
    
        """properties of atoms at each position
        uses array: 
            index 0(3 in the final array with positions): weight, unit: gram
            index 1(4): charge, unit: eV
            index 2(5): epsilon in LJ force with Lithium
            index 3(6): sigma in LJ force with Lithium
                
            Unit Conversion: 1eV= 1.602*10**(-19) J
        """
    
        c_properties_list = np.array([[1.9944*10**(-23),0, 0.105,3.851],[1.9944*10**(-23),0, 0.105,3.851],[1.9944*10**(-23),0, 0.105,3.851],[1.9944*10**(-23),0, 0.105,3.851]])
        fixed_UC_prop = c_properties_list
        UC_ary = np.concatenate((fixed_UC_pos, fixed_UC_prop),axis = 1)
        c_col = np.array([0.1,0.1,0.1,0.1])
        fixed_UC_color = c_col
        adjustedlatticeparam = np.array([self.lattice_params[0]-self.lattice_params[1]*np.cos(np.pi/3),self.lattice_params[1]*np.sin(np.pi/3),self.lattice_params[2]])
        """real space dimension of lattice"""
        dim = np.multiply(adjustedlatticeparam,numofcells)
        """lattice manipulation into cubix"""
        lat , self.col = self.lattice_construction(numofcells, UC_ary, dim, fixed_UC_color)
        lat = self.shift_to_center(lat, dim)
        print(dim)
        lat= self.minimum_image(lat,dim)
        lat = self.removedup(lat)
        print(dim)
        self.mindim = boxdim
        self.latticearray, self.col = self.crop_cubic(lat,self.mindim, self.col)

unitcells = [6,6,3]
boxdimension =7
lat = lattice()
lat.init(unitcells,boxdimension)
a = lat.latticearray
b = lat.col
dim = lat.mindim
print(a.shape)
print(dim)

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
