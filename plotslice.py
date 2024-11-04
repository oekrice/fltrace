#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from fltrace import trace_fieldlines

#matplotlib.rcParams['text.usetex'] = True


if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0


fig_width = 15#1.0*(513.11743/72)


nx = 64
ny = 64
nz = 64

xs = np.linspace(-12,12, nx+1)
ys = np.linspace(-12,12, ny+1)
zs = np.linspace(-24.0/ny,24, nz+1)

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

xc = 0.5*(xs[1:] + xs[:-1])
yc = 0.5*(ys[1:] + ys[:-1])
zc = 0.5*(zs[1:] + zs[:-1])

dx = xs[1] - xs[0]
dy = ys[1] - ys[0]
dz = zs[1] - zs[0]

class Grid():
    def __init__(self):
        self.x0 = xs[0]; self.x1 = xs[-1]
        self.y0 = ys[0]; self.y1 = ys[-1]
        self.z0 = zs[0]; self.z1 = zs[-1]
        self.nx = nx ; self.ny = ny; self.nz = nz

i = 0

remote_flag = 1

if len(sys.argv) > 2:
    i = int(sys.argv[2])

for plot_num in range(1):

    if remote_flag:
        data_directory = './Data/'
    else:
        data_directory = '/home/grads/trcn27/rdata/lare3d_jet/'

    slice_index = ny//2
    i = plot_num
    wait = 0
    fname = '%s%04d.nc' % (data_directory, i)
    print('Making plot', i, 'fname', fname)

    try:
        data = netcdf_file(fname, 'r', mmap=False)
        print('File', fname, 'found')

    except:
        print('File', fname, 'not found')
        continue

    bx = np.zeros((nx+1,ny+2,nz+2))
    by = np.zeros((nx+2,ny+1,nz+2))
    bz = np.zeros((nx+2,ny+2,nz+1))

    bx[:,1:-1,1:-1] = np.swapaxes(data.variables['bx'][:],0,2)
    by[1:-1,:,1:-1] = np.swapaxes(data.variables['by'][:],0,2)
    bz[1:-1,1:-1,:] = np.swapaxes(data.variables['bz'][:],0,2)

    en = np.zeros((nx+2,ny+2,nz+2))
    rho = np.zeros((nx+2,ny+2,nz+2))

    en[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['en'][:],0,2)
    rho[1:-1,1:-1,1:-1] = np.swapaxes(data.variables['rho'][:],0,2)

    vx = np.zeros((nx+1,ny+1,nz+1))
    vy = np.zeros((nx+1,ny+1,nz+1))
    vz = np.zeros((nx+1,ny+1,nz+1))

    vx = np.swapaxes(data.variables['vx'][:],0,2)
    vy = np.swapaxes(data.variables['vy'][:],0,2)
    vz = np.swapaxes(data.variables['vz'][:],0,2)

    jx = np.zeros((nx+2,ny+1,nz+1))
    jy = np.zeros((nx+1,ny+2,nz+1))
    jz = np.zeros((nx+1,ny+1,nz+2))

    jx[1:-1,:,:] = np.swapaxes(data.variables['jx'][:],0,2)
    jy[:,1:-1,:] = np.swapaxes(data.variables['jy'][:],0,2)
    jz[:,:,1:-1] = np.swapaxes(data.variables['jz'][:],0,2)

    pr = rho*en*(2/3)

    data.close()

    def current_test(bx,by,bz):
        jx = (bz[1:-1,1:,:] - bz[1:-1,:-1,:])/dy - (by[1:-1,:,1:] - by[1:-1,:,:-1])/dz
        jy =  (bx[:,1:-1,1:] - bx[:,1:-1,:-1])/dz - (bz[1:,1:-1,:] - bz[:-1,1:-1,:])/dx
        jz =  (by[1:,:,1:-1] - by[:-1,:,1:-1])/dx - (bx[:,1:,1:-1] - bx[:,:-1,1:-1])/dy

        print(np.max(jx), np.max(jy), np.max(jz))

    def magfield(bx, by, bz):
        bx1 = 0.5*(bx[1:,slice_index,1:-1] + bx[:-1,slice_index,1:-1])
        by1 = 0.5*(by[1:-1,slice_index,1:-1] + by[1:-1,slice_index,1:-1])
        bz1 = 0.5*(bz[1:-1,slice_index,1:] + bz[1:-1,slice_index,:-1])
        return 0.5*(bx1**2 + by1**2+ bz1**2)

    if np.max(magfield(bx,by,bz)) > 1e-6:
        beta = 4*np.pi*pr[1:-1,slice_index,1:-1].T/magfield(bx,by,bz).T
    else:
        beta = 0.0*pr[1:-1,slice_index,1:-1].T

    current_test(bx, by, bz)

    if True:
        trace_fieldlines(Grid(),bx,by,bz,save=plot_num,plot_vista = True, plot_notvista = False)

    if False:
        fig, axs = plt.subplots(4,4, figsize = (10,10))

        def find_a(bx, bz):   #find the vector potential a from the (hopefully) divergence-free magnetic fields
            a = np.zeros((nx, nz))
            for j in range(1,nz):
                a[0,j] = a[0,j-1] -dz*bx[0,slice_index,j-1]

            for i in range(1,nx):
                a[i,0] = a[i-1,0] + dx*bz[i-1,slice_index,0]
                for j in range(1,nz):
                    a[i,j] = a[i,j-1] - dz*bx[i,slice_index,j-1]
            return a

        a = find_a(bx, by)

        im = axs[0,0].pcolormesh(xc,zs,bx[1:-1,slice_index,1:-1].T,vmin=-np.max(np.abs(bx[1:-1,slice_index,1:-1])), vmax = np.max(np.abs(bx[1:-1,slice_index,1:-1])), cmap ='seismic')
        plt.colorbar(im, ax=axs[0,0])
        axs[0,0].set_title('Bx')

        im = axs[0,1].pcolormesh(xs,zs,by[1:-1,slice_index,1:-1].T,vmin=-np.max(np.abs(by[1:-1,slice_index,1:-1])), vmax = np.max(np.abs(by[1:-1,slice_index,1:-1])), cmap ='seismic')
        plt.colorbar(im, ax=axs[0,1])
        axs[0,1].set_title('By')

        im = axs[0,2].pcolormesh(xs,zc,bz[1:-1,slice_index,1:-1].T,vmin=-np.max(np.abs(bz[1:-1,slice_index,1:-1])), vmax = np.max(np.abs(bz[1:-1,slice_index,1:-1])), cmap ='seismic')
        plt.colorbar(im, ax=axs[0,2])
        axs[0,2].set_title('Bz')

        im = axs[1,3].pcolormesh(xs,zs,en[1:-1,slice_index,1:-1].T,vmin = 0.0)
        plt.colorbar(im, ax=axs[1,3])
        axs[1,3].set_title('Internal Energy')

        im = axs[1,0].pcolormesh(xs,zs,pr[1:-1,slice_index,1:-1].T,vmin = 0.0)
        plt.colorbar(im, ax=axs[1,0])
        axs[1,0].set_title('Pressure')

        im = axs[1,1].pcolormesh(xs,zs,rho[1:-1,slice_index,1:-1].T,vmin = 0.0)
        plt.colorbar(im, ax=axs[1,1])
        axs[1,1].set_title('Density')

        im = axs[1,2].pcolormesh(xc,zc,vz[1:-1,slice_index,1:-1].T,vmin=-np.max(np.abs(vz)), vmax = np.max(np.abs(vz)),cmap='seismic')
        plt.colorbar(im, ax=axs[1,2])
        axs[1,2].set_title('Vertical Velocity')

        im = axs[0,3].pcolormesh(xc,zc,beta,vmin=0.0, vmax = np.max(np.abs(beta)))
        plt.colorbar(im, ax=axs[0,3])
        axs[0,3].set_title('Plasma Beta')

        im = axs[2,0].pcolormesh(xc,yc,vx[1:-1,1:-1,0])
        plt.colorbar(im, ax=axs[2,0])
        axs[2,0].set_title('Surface Vx')

        im = axs[2,1].pcolormesh(xc,yc,vy[1:-1,1:-1,0])
        plt.colorbar(im, ax=axs[2,1])
        axs[2,1].set_title('Surface Vy')

        im = axs[2,2].pcolormesh(xc,yc,bz[1:-1,1:-1,0])
        plt.colorbar(im, ax=axs[2,2])
        axs[2,2].set_title('Surface Bz')

        print(np.shape(jx), np.shape(xs), np.shape(zc))
        im = axs[3,0].pcolormesh(xs,zc,jx[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[3,0])
        axs[3,0].set_title('Current jx')

        im = axs[3,1].pcolormesh(xc,zc,jy[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[3,1])
        axs[3,1].set_title('Current jy')

        im = axs[3,2].pcolormesh(xc,zs,jz[1:-1,slice_index,1:-1].T)
        plt.colorbar(im, ax=axs[3,2])
        axs[3,2].set_title('Current jz')


        print(beta[slice_index, slice_index])
        plt.tight_layout()
        #plt.show()
        plt.savefig('plots/a%04d' % plot_num)
        plt.close()

    if False:
        os.system('ffmpeg -y -framerate 20 -i ./plots/a%04d.png -b:v 10M  sideways2.mp4')
