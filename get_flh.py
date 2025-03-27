#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np 
from matplotlib import pyplot as plt
import field_line_topology as flt
from streamtracer import StreamTracer, VectorGrid
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import simpson
#import waveletRoutines as wr

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import netcdf_file
import os
import sys
# Script called from the field line tracer to get a distribution of the FLH on the lower surface.

class FLH():
    def __init__(self, grid, closed_boundaries = False):

        high_resolution = False
        #Get grid and stuff in here so don't have to do it twice
        bxOg = grid.bx
        byOg = grid.by
        bzOg = grid.bz

        bx = 0.5*(bxOg[1:,:,:] + bxOg[:-1,:,:])
        by = 0.5*(byOg[:,1:,:] + byOg[:,:-1,:])
        bz = 0.5*(bzOg[:,:,1:] + bzOg[:,:,:-1])

        bz = bz - np.sum(bz[:,:,0])/np.size(bz[:,:,0])

        xv = grid.xc
        yv = grid.yc
        zv = grid.zc

        X, Y = np.meshgrid(xv, yv, indexing='ij')
        # Flatten the grid arrays to form the input to the interpolator
        points = np.vstack([X.ravel(), Y.ravel()]).T
        dx = xv[1]-xv[0]
        dy = yv[1]-yv[0]
        dz = zv[1]-zv[0]
        dA = dx*dy
        grid_spacing = [dx,dy,dz]
        grid_ncells = [grid.nx, grid.ny, grid.nz]

        bField = flt.createSingleField(bx,by,bz)

        AField = flt.getAFastSingle(bField,grid_ncells,grid_spacing)

        if not closed_boundaries:
            AField = flt.adjust_A_boundary_flux(AField, bField, xv, yv, zv)

        #Account for nonzero net flux through the lower boundary.
        bzConst = np.sum(bz[:,:,0])/(dA*(grid.nx)*(grid.ny))
        AConst = flt.AConst(bzConst,points,dA,[grid.nx, grid.ny, grid.nz])
        AField = AField + AConst
    
            
        bField_test = flt.curl(AField,grid_spacing)

        usf = flt.unitSpeedField(bField.copy(),0.5)  #transforms bz to be 'unit speed' in the z direction
        fig, axs = plt.subplots(1,3)
        ax = axs[0]
        ax.pcolormesh(usf[:,:,0,0])
        ax.set_aspect('equal')
        ax.set_title('ux')
        ax = axs[1]
        ax.pcolormesh(usf[:,:,0,1])
        ax.set_aspect('equal')
        ax.set_title('uy')
        ax = axs[2]
        ax.pcolormesh(usf[:,:,0,2])
        ax.set_aspect('equal')
        ax.set_title('uz')
        plt.tight_layout()
        plt.close()
        BUnit = flt.addDivergenceCleaningTerm(usf,grid_ncells,grid_spacing)   #Returns unit speed field which is roughly divergence-free
        divTestA = flt.divergence(usf,grid_spacing)
        divTestB = flt.divergence(BUnit,grid_spacing)
        fig, axs = plt.subplots(1,2)
        ax = axs[0]
        im = ax.pcolormesh(divTestA[:,:,0],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('ux')
        ax = axs[1]
        im = ax.pcolormesh(divTestB[:,:,0],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('uy')
        plt.close()
        AWind = flt.getAFastSingle(BUnit,grid_ncells,grid_spacing)  #Returns A for the winding


        #Account for nonzero net flux through the lower boundary.
        bzUnitConst = np.sum(BUnit[:,:,0,2])/(dA*(grid.nx)*(grid.ny))
        AUnitConst = flt.AConst(bzUnitConst,points,dA,[grid.nx, grid.ny, grid.nz])
        AWind = AWind + AUnitConst
        
        if not closed_boundaries:
            AWind = flt.adjust_A_boundary_flux(AWind, BUnit, xv, yv, zv)
            
        aWindTest = flt.curl(AWind,grid_spacing)
        fig, axs = plt.subplots(2,2)
        ax = axs[0][0]
        im = ax.pcolormesh(BUnit[:,:,0,0],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('ux')
        ax = axs[0][1]
        im = ax.pcolormesh(aWindTest[:,:,0,0],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('uy')
        ax = axs[1][0]
        im = ax.pcolormesh(BUnit[:,:,0,1],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('ux')
        ax = axs[1][1]
        im =ax.pcolormesh(aWindTest[:,:,0,1],vmin=-5,vmax=5)
        fig.colorbar(im,ax=ax)
        ax.set_aspect('equal')
        ax.set_title('uy')
        plt.close()
            
        curlField = flt.curl(bField,grid_spacing)

        self.flh_density = flt.getFLHDenSingle(bField,AField,signed = True)
        self.winding_density = flt.getFLHDenSingle(usf,AWind,signed = True)
        self.twist_density = flt.twistDen(bField,curlField,0.001)









