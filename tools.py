#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for fortran field line tracer

Uses FLH density on the photosphere to determine which lines to trace.

Needs routine from emergence_comparisons to do that. Bit of a mess...
"""
import numpy as np
from scipy.io import netcdf_file

def interpolate_to_staggered_grid(bx_in, by_in, bz_in, swapaxes = False):
    #Assuming that the dimensions of the input field are all the same, interpolates to an interior staggered grid
    #Check grids are all fine
    if swapaxes:
        bx_in = np.swapaxes(bx_in, 0, 2)
        by_in = np.swapaxes(by_in, 0, 2)
        bz_in = np.swapaxes(bz_in, 0, 2)

    nx = bx_in.shape[0]; ny = bx_in.shape[1]; nz = bx_in.shape[2]
    if not(by_in.shape[0] == nx and by_in.shape[1] == ny and by_in.shape[2] == nz):
        raise Exception("Input files don't have the same dimensions")
    if not(bz_in.shape[0] == nx and bz_in.shape[1] == ny and bz_in.shape[2] == nz):
        raise Exception("Input files don't have the same dimensions")

    bx_out = np.zeros((bx_in.shape[0]+1, bx_in.shape[1], bx_in.shape[2]))
    by_out = np.zeros((by_in.shape[0], by_in.shape[1] + 1, by_in.shape[2]))
    bz_out = np.zeros((bz_in.shape[0], bz_in.shape[1], bz_in.shape[2] + 1))

    bx_out[1:-1,:,:] = 0.5*(bx_in[1:,:,:] + bx_in[:-1,:,:])
    by_out[:,1:-1,:] = 0.5*(by_in[:,1:,:] + by_in[:,:-1,:])
    bz_out[:,:,1:-1] = 0.5*(bz_in[:,:,1:] + bz_in[:,:,:-1])

    bx_out[0,:,:] = bx_out[1,:,:]; bx_out[-1,:,:] = bx_out[-2,:,:]
    by_out[:,0,:] = by_out[:,1,:]; by_out[:,-1,:] = by_out[:,-2,:]
    bz_out[:,:,0] = bz_out[:,:,1]; bz_out[:,:,-1] = bz_out[:,:,-2]
    return bx_out, by_out, bz_out

def check_staggered_grid(bx_in, by_in, bz_in):
    #One dimension in each of the coordinates should be one cell longer than the others.
    #OR theoretically could be one shorter.
    if bx_in.shape[0] == bx_in.shape[1]:
        xcoord = 2
    elif bx_in.shape[1] == bx_in.shape[2]:
        xcoord = 0
    else:
        xcoord = 1
    if by_in.shape[0] == by_in.shape[1]:
        ycoord = 2
    elif by_in.shape[1] == by_in.shape[2]:
        ycoord = 0
    else:
        ycoord = 1
    #Determine if axes need rearranging -- hopefully just a straight swap...
    if xcoord == 0 and ycoord == 1:
        pass
    elif xcoord == 0 and ycoord == 2:
        bx_in = np.swapaxes(bx_in, 1, 2)
        by_in = np.swapaxes(by_in, 1, 2)
        bz_in = np.swapaxes(bz_in, 1, 2)
    elif xcoord == 1 and ycoord == 0:
        bx_in = np.swapaxes(np.swapaxes(bx_in, 1, 2), 0, 2)
        by_in = np.swapaxes(np.swapaxes(by_in, 1, 2), 0, 2)
        bz_in = np.swapaxes(np.swapaxes(bz_in, 1, 2), 0, 2)
    elif xcoord == 1 and ycoord == 2:
        bx_in = np.swapaxes(bx_in, 0, 1)
        by_in = np.swapaxes(by_in, 0, 1)
        bz_in = np.swapaxes(bz_in, 0, 1)
    elif xcoord == 2 and ycoord == 1:
        bx_in = np.swapaxes(bx_in, 0, 2)
        by_in = np.swapaxes(by_in, 0, 2)
        bz_in = np.swapaxes(bz_in, 0, 2)
    elif xcoord == 2 and ycoord == 0:
        bx_in = np.swapaxes(np.swapaxes(bx_in, 0, 1), 0,2)
        by_in = np.swapaxes(np.swapaxes(by_in, 0, 1), 0,2)
        bz_in = np.swapaxes(np.swapaxes(bz_in, 0, 1), 0,2)
    else:
        raise Exception("Can't make sense of input dimensions...")

    #Check if need to remove extra cells
    if bx_in.shape[1] > bx_in.shape[0]:
        bx_in = bx_in[:,1:-1,1:-1]
        by_in = by_in[1:-1,:,1:-1]
        bz_in = bz_in[1:-1,1:-1,:]

    return bx_in, by_in, bz_in

def make_coordinate_arrays(Grid):
    #If coordinate arrays don't exist, make them. If they do, check they make sense.
    #Replace things in place for neatness purposes
    if Grid.xs is None:  #Establish a grid with equal cubular cells, and zero as the base z coordinate
        Grid.x0 = -1; Grid.x1 = 1
        Grid.y0 = Grid.x0*Grid.ny/Grid.nx; Grid.y1 = Grid.x1*Grid.ny/Grid.nx
        Grid.z0 = 0.0; Grid.z1 = 2.0*Grid.x1*Grid.nz/Grid.nx

        Grid.xs = np.linspace(Grid.x0,Grid.x1,Grid.nx+1)
        Grid.ys = np.linspace(Grid.y0,Grid.y1,Grid.ny+1)
        Grid.zs = np.linspace(Grid.z0,Grid.z1,Grid.nz+1)

        Grid.xc = 0.5*(Grid.xs[1:] + Grid.xs[:-1])
        Grid.yc = 0.5*(Grid.ys[1:] + Grid.ys[:-1])
        Grid.zc = 0.5*(Grid.zs[1:] + Grid.zs[:-1])

        Grid.dx = Grid.xs[1] - Grid.xs[0]
        Grid.dy = Grid.ys[1] - Grid.ys[0]
        Grid.dz = Grid.zs[1] - Grid.zs[0]

    else:
        #Check things are the right size. If imported these will probably actually be the grid centres
        if (len(Grid.xs) == Grid.nx and len(Grid.ys) == Grid.ny and len(Grid.zs) == Grid.nz):
            #This is on the centres
            Grid.xc = Grid.xs
            Grid.yc = Grid.ys
            Grid.zc = Grid.zs
            Grid.dx = Grid.xc[1] - Grid.xc[0]; Grid.dy = Grid.yc[1] - Grid.yc[0]; Grid.dz = Grid.zc[1] - Grid.zc[0]
            Grid.x0 = Grid.xc[0] - Grid.dx/2; Grid.x1 = Grid.xc[-1] + Grid.dx/2
            Grid.y0 = Grid.yc[0] - Grid.dy/2; Grid.y1 = Grid.yc[-1] + Grid.dy/2
            Grid.z0 = Grid.zc[0] - Grid.dz/2; Grid.z1 = Grid.zc[-1] + Grid.dz/2
            Grid.xs = np.linspace(Grid.x0,Grid.x1,Grid.nx+1)
            Grid.ys = np.linspace(Grid.y0,Grid.y1,Grid.ny+1)
            Grid.zs = np.linspace(Grid.z0,Grid.z1,Grid.nz+1)

        elif (len(Grid.xs) == Grid.nx + 1 and len(Grid.ys) == Grid.ny + 1 and len(Grid.zs) == Grid.nz + 1):
            #This is on the points
            Grid.xc = 0.5*(Grid.xs[1:] + Grid.xs[:-1])
            Grid.yc = 0.5*(Grid.ys[1:] + Grid.ys[:-1])
            Grid.zc = 0.5*(Grid.zs[1:] + Grid.zs[:-1])
            Grid.dx = Grid.xs[1] - Grid.xs[0]
            Grid.dy = Grid.ys[1] - Grid.ys[0]
            Grid.dz = Grid.zs[1] - Grid.zs[0]
            Grid.x0 = Grid.xs[0]; Grid.x1 = Grid.xs[-1]
            Grid.y0 = Grid.ys[0]; Grid.y1 = Grid.ys[-1]
            Grid.z0 = Grid.zs[0]; Grid.z1 = Grid.zs[-1]
        else:
            raise Exception("Imported coordinates don't match array sizes")

    if True:
        #Shift such that the lower z coordinate is zero
        Grid.zc = Grid.zc - Grid.z0
        Grid.zs = Grid.zs - Grid.z0
        Grid.z1 = Grid.z1 - Grid.z0
        Grid.z0 = 0
    return

def calculate_current(Grid):
    Grid.jx = np.zeros((Grid.nx,Grid.ny+1,Grid.nz+1))
    Grid.jy = np.zeros((Grid.nx+1,Grid.ny,Grid.nz+1))
    Grid.jz = np.zeros((Grid.nx+1,Grid.ny+1,Grid.nz))

    Grid.jx[:,1:-1,1:-1] =  (Grid.bz[:,1:,1:-1] - Grid.bz[:,:-1,1:-1])/Grid.dy - (Grid.by[:,1:-1,1:] - Grid.by[:,1:-1,:-1])/Grid.dz

    Grid.jy[1:-1,:,1:-1] =  (Grid.bx[1:-1,:,1:] - Grid.bx[1:-1,:,:-1])/Grid.dz - (Grid.bz[1:,:,1:-1] - Grid.bz[:-1,:,1:-1])/Grid.dx
    Grid.jz[1:-1,1:-1,:] =  (Grid.by[1:,1:-1,:] - Grid.by[:-1,1:-1,:])/Grid.dx - (Grid.bx[1:-1,1:,:] - Grid.bx[1:-1,:-1,:])/Grid.dy

    return

def save_for_fortran(Grid):
    #Save out to temporary file so Fortran can read it after it has been bodged
    fid = netcdf_file('./tmp/%05d.nc' % (Grid.id), 'a')
    fid.createDimension('xs', Grid.nx+1)
    fid.createDimension('ys', Grid.ny+1)
    fid.createDimension('zs', Grid.nz+1)
    fid.createDimension('xc', Grid.nx)
    fid.createDimension('yc', Grid.ny)
    fid.createDimension('zc', Grid.nz)

    vid = fid.createVariable('xs', 'd', ('xs',))
    vid[:] = Grid.xs
    vid = fid.createVariable('ys', 'd', ('ys',))
    vid[:] = Grid.ys
    vid = fid.createVariable('zs', 'd', ('zs',))
    vid[:] = Grid.zs

    vid = fid.createVariable('xc', 'd', ('xc',))
    vid[:] = Grid.xc
    vid = fid.createVariable('yc', 'd', ('yc',))
    vid[:] = Grid.yc
    vid = fid.createVariable('zc', 'd', ('zc',))
    vid[:] = Grid.zc

    #Transposes are necessary as it's easier to flip here than in Fortran
    vid = fid.createVariable('bx', 'd', ('zc','yc','xs'))
    vid[:] = np.swapaxes(Grid.bx, 0, 2)
    vid = fid.createVariable('by', 'd', ('zc','ys','xc'))
    vid[:] = np.swapaxes(Grid.by, 0, 2)
    vid = fid.createVariable('bz', 'd', ('zs','yc','xc'))
    vid[:] = np.swapaxes(Grid.bz, 0, 2)

    vid = fid.createVariable('jx', 'd', ('zs','ys','xc'))
    vid[:] = np.swapaxes(Grid.jx, 0, 2)
    vid = fid.createVariable('jy', 'd', ('zs','yc','xs'))
    vid[:] = np.swapaxes(Grid.jy, 0, 2)
    vid = fid.createVariable('jz', 'd', ('zc','ys','xs'))
    vid[:] = np.swapaxes(Grid.jz, 0, 2)

    #Now the flh quantities, if necessary. Check they exist first
    if Grid.flh_density is not None:
        vid = fid.createVariable('flh_density', 'd', ('zc','yc','xc'))
        vid[:] = np.swapaxes(Grid.flh_density, 0, 2)
    if Grid.winding_density is not None:
        vid = fid.createVariable('winding_density', 'd', ('zc','yc','xc'))
        vid[:] = np.swapaxes(Grid.winding_density, 0, 2)
    if Grid.twist_density is not None:
        vid = fid.createVariable('twist_density', 'd', ('zc','yc','xc'))
        vid[:] = np.swapaxes(Grid.twist_density, 0, 2)

    fid.close()

    return









