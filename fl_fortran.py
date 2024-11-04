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

#Establish grid parameters (can be read in from elsewhere of course)
run = 0
snap = 0
print_flag = 1

nx = 64
ny = 64
nz = 64

x0 = -12.; x1 = 12.
y0 = -12.; y1 = 12.
z0 = -24.0/ny; z1 = 24

#Establish start points for the field line plotting
max_line_length = 10000
ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
weakness_limit = 1e-3   #Minimum field strength to stop plotting

starts = []
#Trace from the top
nrs = 10; nthetas = 2
ris = np.linspace(3.0/nrs,x0-3.0/nrs,nrs)
tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
for i in range(nrs):
    for j in range(nthetas):
        starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),z1-1e-6])
#And from the bottom (for the interior ones only)
nrs = 20; nthetas = 2
ris = np.linspace(0.5*x0/nrs,x0-0.5*x0/nrs,nrs)
tjs = np.linspace(0+1e-6,2*np.pi*(1-1/nthetas),nthetas)
for i in range(nrs):
    for j in range(nthetas):
        starts.append([ris[i]*np.cos(tjs[j]),ris[i]*np.sin(tjs[j]),1e-6])

nstarts = len(starts)
starts = np.array(starts).reshape(nstarts*3)

variables = np.zeros((30))

variables[0] = run
variables[1] = nx
variables[2] = ny
variables[3] = nz
variables[4] = x0
variables[5] = x1
variables[6] = y0
variables[7] = y1
variables[8] = z0
variables[9] = z1
variables[10] = snap
variables[11] = nstarts
variables[12] = print_flag
variables[13] = max_line_length
variables[14] = ds
variables[15] = weakness_limit

np.savetxt('flparameters.txt', variables)   #variables numbered based on run number (up to 1000)
np.savetxt('starts.txt', starts)   #Coordinates of the start points of each field line (do this in python)

os.system('make')
os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace')
