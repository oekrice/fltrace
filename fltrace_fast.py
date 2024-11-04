#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for fortran field line tracer

"""


import os
import shutil
import numpy as np
import sys
from numpy import random
import time
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
import pyvista as pv
pv.start_xvfb()

class trace_fieldlines():
    def __init__(self):
        #Establish grid parameters (can be read in from elsewhere of course)
        self.run = 0
        self.snap = 0
        self.print_flag = 1
        self.save_number = self.snap
        self.data_root = './Data/'

        self.bz = np.swapaxes(netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False).variables['bz'][:],0,2)
        #Import bz as a test of the resolutions (and for the pyvista plot)
        self.nx = np.shape(self.bz)[0]
        self.ny = np.shape(self.bz)[1]
        self.nz = np.shape(self.bz)[2] - 1

        self.x0 = -12.; self.x1 = 12.
        self.y0 = -12.; self.y1 = 12.
        self.z0 = -24.0/self.nz; self.z1 = 24

        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)

        #Establish start points for the field line plotting
        self.max_line_length = 10000
        self.ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
        self.weakness_limit = 1e-3   #Minimum field strength to stop plotting
        self.line_plot_length = 50  #To save time while plotting, reduce the length of the plotted lines

        #Folder admin
        if not os.path.exists('./fl_data/'):
            os.mkdir('fl_data')
        os.system('rm ./fl_data/flines.nc')
        #Find start points
        self.set_starts()
        #Create runtime variables for fortran
        self.setup_tracer()
        #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
        self.trace_lines_fortran()
        #Plot the field lines (using pyvista)
        if True:
            if not os.path.exists('./plots/'):
                os.mkdir('plots')
            self.plot_vista()

    def plot_vista(self):
        print('Plotting...')
        x, y = np.meshgrid(self.xs, self.ys)
        z = 0*x*y
        surface = pv.StructuredGrid(x, y, z)
        p = pv.Plotter(off_screen=True)
        p.background_color = "black"
        p.add_mesh(surface, scalars= self.bz[:,:,0], show_edges=True,cmap = 'plasma')


        for li, line in enumerate(self.lines):
            line = np.array(line)
            line_length = len(line[line[:,2]<1e6])
            #Thin out the lines (if required)
            if line_length > 0:
                thin_fact = max(int(line_length/self.line_plot_length), 1)
                thinned_line = line[:line_length:thin_fact].copy()
                thinned_line[-1] = line[line_length-1].copy()
            else:
                continue

            line = np.array(thinned_line).tolist()
            doplot = True
            if line_length == 0:
                doplot = False
            elif line[0][2] < 1.0 and line[-1][2] > 20.0:
                doplot = False

            if doplot:

                p.add_mesh(pv.Spline(line, len(line)),color='white')

        p.camera.position = (20.0,40,20.0)
        p.camera.focal_point = (0,0,4)
        p.show(screenshot='plots/b%04d.png' % self.save_number, window_size = (1000,1000))
        print('Plot saved to file plots/b%04d.png' % self.save_number)

    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.
        self.starts = []
        #Trace from the top
        nrs = 10; nthetas = 2
        ris = np.linspace(3.0/nrs,self.x0-3.0/nrs,nrs)
        tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
        for i in range(nrs):
            for j in range(nthetas):
                self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),self.z1-1e-6])
        #And from the bottom (for the interior ones only)
        nrs = 20; nthetas = 10
        ris = np.linspace(0.5*self.x0/nrs,self.x0-0.5*self.x0/nrs,nrs)
        tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
        for i in range(nrs):
            for j in range(nthetas):
                self.starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),1e-6])

        self.nstarts = len(self.starts)
        self.starts = np.array(self.starts).reshape(self.nstarts*3)

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.05 #Tracing 'timestep' as a proportion of the grid size
        weakness_limit = 1e-3   #Minimum field strength to stop plotting
        print_flag = 1  #Print some things as the tracing happens

        variables = np.zeros((30))

        variables[0] = self.run
        variables[1] = self.nx
        variables[2] = self.ny
        variables[3] = self.nz
        variables[4] = self.x0
        variables[5] = self.x1
        variables[6] = self.y0
        variables[7] = self.y1
        variables[8] = self.z0
        variables[9] = self.z1
        variables[10] = self.snap
        variables[11] = self.nstarts
        variables[12] = self.print_flag
        variables[13] = self.max_line_length
        variables[14] = self.ds
        variables[15] = self.weakness_limit

        np.savetxt('./fl_data/flparameters.txt', variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('./fl_data/starts.txt', self.starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):
        os.system('make')
        if os.uname()[1] == 'brillouin.dur.ac.uk':
            os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace')
        elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
            os.system('mpiexec -np 1 ./bin/fltrace')
        else:
            os.system('mpirun -np 1 ./bin/fltrace')

        try:
            data = netcdf_file('./fl_data/flines.nc', 'r', mmap=False)
            print('Field lines found')

        except:
            print('File not found')

        self.lines = np.swapaxes(data.variables['lines'][:],0,2)

trace_fieldlines()
