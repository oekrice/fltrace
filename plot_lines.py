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
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter
from tools import interpolate_to_staggered_grid, check_staggered_grid, make_coordinate_arrays, calculate_current, save_for_fortran
from get_flh import FLH
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator

import matplotlib
matplotlib.rcParams['text.usetex'] = True

'''
Tools for 3D magnetic field analysis.
Requires input as a .nc file, either on a staggered grid or not.
If there is dimension data it will use that, and if not will make something up.
Calculates FLH, winding and twist in python
Saved modified magnetic field to a temporary file
Does field line tracing in fortran, which calculated the quantities integrated along field lines, and white light emissions.
Reads back into python and plots in 'plot_emissions'
IMPORTANT -- WILL NEED TO CHANGE THE STUFF AT THE TOP OF THE MAKEFILE DEPENDING ON THE MACHINE
Will save plots to ./plots/ folder as pngs.

This one is for using existing emissions to decide which field lines to plot, then output as a 3D plot
'''

for id in range(0,1):   #Loop for multiple runs

    id = 61200   #This should be the same ID as used in fltrace for the emissions:
    #eg. the output is fl_data/emiss61200.nc if id = 61200.
    #LEADING ZEROS ARE NECESSARY

    #Define parameters for plotting:
    input_fname = 'Bout__2021.0621.061200.nc'

    nlines = 5000                            #Approx number of field lines to trace
    show = True                            #Brings up pyvista interactively. If false will save to /plots. Can't do both for reasons I don't understand.
    justplot = False                       #If true, finds existing data and just plots it. If you want to just tweak the plots without running everything again.
    remove_emission_files = False               #Removes the emission files after plotting
    swapaxes = True                             #Swaps x and z axes of the imported file

    plot_angle = 0.0   #Angle in radians from the default for the final screenshot
    class Fltrace():
        def __init__(self, input_fname, id = 0, nlines = 10000, show = True, do_flh = False, swapaxes = False, remove_files = False):

            #Some parameters
            self.max_line_length = 10000
            self.ds = 0.05 #Tracing 'timestep' as a proportion of the self size
            self.weakness_limit = 1e-5  #Minimum field strength to stop plotting
            self.line_plot_length = 100  #To save time while plotting, reduce the length of the plotted lines
            self.nlines = nlines/2

            self.input_fname = input_fname
            self.id = id
            self.remove_files = remove_files
            #In case do_flh isn't on, this stops an error happening later
            self.flh_density = None
            self.winding_density = None
            self.twist_density = None

            #Machine-pecific things for saving out etc.
            data_directory = ('./input/')

            #Create some folders if necessary
            if not os.path.isdir(data_directory):
                os.mkdir(data_directory)
            if not os.path.isdir('./plots/'):
                os.mkdir('./plots/')
            if not os.path.isdir('./fl_data/'):
                os.mkdir('./fl_data/')
            if not os.path.isdir('./tmp/'):
                os.mkdir('./tmp/')

            self.show = show
            #Establish self parameters (can be read in from elsewhere of course)
            self.print_flag = show

            self.option = 2   #tracing a plotting options (1 for jet, 2 for emergence)

            self.data_source = 0.
            self.data_root = data_directory

            #Establish random seeds for field line plotting
            if not os.path.isfile('./fl_data/start_seeds.txt'):
                self.start_seeds = random.rand(1000**2)
                np.savetxt('./fl_data/start_seeds.txt', self.start_seeds, delimiter = ',')
            else:
                self.start_seeds = np.loadtxt('./fl_data/start_seeds.txt', delimiter = ',')
            data = netcdf_file('%s%s' % (self.data_root, input_fname), 'r', mmap=False)

            print('Using data from file:', '%s%s' % (self.data_root, input_fname))
            #Copy file so fotran can read it. OR, write it out proper like.
            os.system('scp -r %s%s ./tmp/%05d.nc' % (self.data_root, input_fname, self.id))
            self.xs = None; self.ys = None; self.zs = None
            #Establish if on a staggered self or not
            bx_input = data.variables['bx'][:]
            by_input = data.variables['by'][:]
            bz_input = data.variables['bz'][:]
            print('Import dimensions', bx_input.shape, by_input.shape, bz_input.shape)

            if bx_input.shape == by_input.shape:
                #This is not on a staggered grid -- interpolate to one for tracing?
                #Assume it is on centres when imported.
                #Will assume axes don't need flipping unless it's otherwise obvious
                self.bx, self.by, self.bz = interpolate_to_staggered_grid(bx_input, by_input, bz_input, swapaxes = swapaxes)
            else:
                #Just check the coordinates make sense, and flip axes if necessary
                self.bx, self.by, self.bz = check_staggered_grid(bx_input, by_input, bz_input)

            #Also need b on grid centres to do the flh calculations. So do that.
            print('Adjusted dimensions', self.bx.shape, self.by.shape, self.bz.shape)

            #Get number of grid cells
            self.nx = np.shape(self.bx)[0] - 1
            self.ny = np.shape(self.bx)[1]
            self.nz = np.shape(self.bx)[2]
            try:
                self.xs = data.variables['x'][:]
                self.ys = data.variables['y'][:]
                self.zs = data.variables['z'][:]
            except:
                pass
            data.close()

            make_coordinate_arrays(self)

            if False: #Plot slice of each coordinate to check this looks right. If it doesn't, try swapping the import axes?
                slice_number = 0
                fig, axs = plt.subplots(1,3)
                ax = axs[0]
                ax.pcolormesh(self.bx[:,:,slice_number])
                ax.set_aspect('equal')
                ax.set_title('Bx')
                ax = axs[1]
                ax.pcolormesh(self.by[:,:,slice_number])
                ax.set_aspect('equal')
                ax.set_title('By')
                ax = axs[2]
                ax.pcolormesh(self.bz[:,:,slice_number])
                ax.set_aspect('equal')
                ax.set_title('Bz')
                plt.suptitle('Imported magnetic field at height %.1f' % self.zs[slice_number])
                plt.tight_layout()
                plt.show()

            print('Coordinate ncells', len(self.xc), len(self.yc), len(self.zc))
            print('Coordinate boundaries', [self.x0, self.x1], [self.y0, self.y1], [self.z0, self.z1])

            calculate_current(self)   #This may have been imported, but also might not have been so just calculate it anyways

            print('Saving for use by Fortran')
            save_for_fortran(self)
            print('Saved to temporary file')

            #Folder admin
            if not os.path.exists('./fl_data/'):
                os.mkdir('fl_data')

            #Find start points. Assumed to be on a regularly-spaced grid
            self.set_starts()

        def trace_lines(self):

            #Create runtime variables for fortran and save out to tmp folder
            self.setup_tracer()
            #Do the tracing
            self.trace_lines_fortran()

            #Remove temporary files
            if False:
                os.system('rm ./fl_data/flparameters%05d.txt' % self.id)
                os.system('rm ./fl_data/starts%05d.txt' % self.id)
                os.system('rm ./tmp/%05d.nc' % (self.id))

        def interpolate_surface_array(self, array_in):
            #For when the surface you want is not the correct resolution. Does 2dinterp or something.
            #Target resolution self.nx, self.ny

            xs_import = np.linspace(0, 1, array_in.shape[0])
            ys_import = np.linspace(0, 1, array_in.shape[1])

            xs_out = np.linspace(0,1,self.nx)
            ys_out = np.linspace(0,1,self.ny)

            X, Y = np.meshgrid(xs_out, ys_out, indexing = 'ij')
            fn = RegularGridInterpolator((xs_import, ys_import), array_in, bounds_error = False, method = 'linear', fill_value = 0.0)
            array_out = fn((X,Y))   #Difference now interpolated to the new grid

            return array_out



        def set_starts(self):
            #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.

            data = netcdf_file('./fl_data/emiss%05d.nc' % self.id, 'r', mmap=False)

            self.flh_array = data.variables['flh_array'][:]
            self.winding_array = data.variables['winding_array'][:]
            self.twist_array = data.variables['twist_array'][:]
            self.mag_array = data.variables['surface_array'][:]
            self.current_array = data.variables['current_array'][:]

            data.close()

            nlines = 500000   #Number from the original integration, not this one
            nx_lines = int(np.sqrt(nlines)*(self.x1-self.x0)/(self.y1 - self.y0))
            ny_lines = int(nlines/nx_lines)

            flh_mesh = np.zeros((nx_lines, ny_lines))
            winding_mesh = np.zeros((nx_lines, ny_lines))
            twist_mesh = np.zeros((nx_lines, ny_lines))
            mag_mesh = np.zeros((nx_lines, ny_lines))
            current_mesh = np.zeros((nx_lines, ny_lines))

            count = 0
            for i in range(nx_lines):
                for j in range(ny_lines):
                    flh_mesh[i,j] = self.flh_array[count]
                    winding_mesh[i,j] = self.winding_array[count]
                    twist_mesh[i,j] = self.twist_array[count]
                    mag_mesh[i,j] = self.mag_array[count]
                    current_mesh[i,j] = self.current_array[count]
                    count += 1

            #_______________________________________________________________________________________________________________
            reference_array = np.sign(mag_mesh)*winding_mesh   #CHANGE THIS IF YOU WANT SOMETHING DIFFERENT
            #reference_array = flh_mesh
            #reference_array = mag_mesh

            #reference_array = np.ones(winding_mesh.shape)



            self.surface_array = self.interpolate_surface_array(reference_array)   #Distribution of surface somethingorother
            #Want to plot the lines based on where the flh differences are highest

            max_surface = np.max(np.abs(self.surface_array + 1)) + 1e-6

            nlines = self.nlines

            alpha = 2.0
            alphasum = np.sum(np.abs(self.surface_array + 1)**alpha)
            pb = max_surface**alpha*nlines/alphasum

            self.starts = []

            #Add capcbility for higher-resolution base
            xscale = np.shape(self.surface_array)[0]/self.nx
            yscale = np.shape(self.surface_array)[1]/self.ny

            nxl = np.shape(self.surface_array)[0]
            nyl = np.shape(self.surface_array)[1]

            self.xsl = np.linspace(self.xs[0], self.xs[-1], nxl)
            self.ysl = np.linspace(self.ys[0], self.ys[-1], nyl)

            xcl = 0.5*(self.xsl[1:] + self.xsl[:-1])
            ycl = 0.5*(self.ysl[1:] + self.ysl[:-1])

            cellcount = 0
            for i in range(nxl-1):  #run through lower surface cells
                for j in range(nyl-1):
                    prop = np.abs(self.surface_array[i,j] + 1)/max_surface
                    #if prop > 0.9:
                    if self.start_seeds[cellcount] < pb*prop**alpha:
                        self.starts.append([xcl[i],ycl[j],0.0])
                    cellcount += 1

            print('Tracing', len(self.starts), 'lines')

            self.nstarts = len(self.starts)
            self.starts = np.array(self.starts).reshape(self.nstarts*3)

        def setup_tracer(self):
            #Output runtime variables to be read-in to Fortran code
            max_line_length = 10000
            ds = 0.1 #Tracing 'timestep' as a proportion of the self size
            print_flag = 1  #Print some things as the tracing happens
            save_all = 1

            self.nx_out = 1024   #Resolutions of the output data
            self.ny_out = int(self.nx_out*self.ny/self.nx)
            self.nz_out = int(self.nx_out*self.nz/self.nx)

            variables = np.zeros((30))

            variables[0] = 0
            variables[1] = self.nx
            variables[2] = self.ny
            variables[3] = self.nz
            variables[4] = self.x0
            variables[5] = self.x1
            variables[6] = self.y0
            variables[7] = self.y1
            variables[8] = self.z0
            variables[9] = self.z1
            variables[10] = self.id
            variables[11] = self.nstarts
            variables[12] = self.print_flag
            variables[13] = self.max_line_length
            variables[14] = self.ds
            variables[15] = self.weakness_limit
            variables[16] = self.data_source
            variables[17] = 0
            variables[18] = self.nx_out
            variables[19] = self.ny_out
            variables[20] = self.nz_out
            variables[21] = print_flag
            variables[22] = save_all


            np.savetxt('./fl_data/flparameters%05d.txt' % self.id, variables)   #variables numbered based on run number (up to 1000)
            np.savetxt('./fl_data/starts%05d.txt' % self.id, self.starts)   #Coordinates of the start points of each field line (do this in python)

        def trace_lines_fortran(self):
            os.system('make')
            os.system('./bin/fltrace %d' % self.id)

        def plot_fieldlines(self):
            #Plots field lines in pyvista
            data = netcdf_file('./fl_data/flines%05d.nc' % (self.id), 'r', mmap=False)

            self.lines = np.swapaxes(data.variables['lines'][:],0,2)

            data.close()

            if show:
                p = pv.Plotter(off_screen=False, shape = (1,1))
            else:
                p = pv.Plotter(off_screen=True, shape = (1,1))

            for li, line in enumerate(self.lines):
                line_length = 0
                for k in range(len(line)):
                    if line[k,2] < 1e6:
                        line_length += 1
                    else:
                        break
                line = line[:line_length,:]
                #Thin out the lines (if required)
                if line_length > 1:
                    thin_fact = max(int(line_length/self.line_plot_length), 1)
                    thinned_line = line[:line_length:thin_fact].copy()
                    thinned_line[-1] = line[line_length-1].copy()
                else:
                    continue

                line = np.array(thinned_line).tolist()
                doplot = True

                if doplot:
                    p.add_mesh(pv.Spline(line, len(line)),color='white',line_width=0.1)

            z_photo = int((self.nz)*(10.0 - self.z0)/(self.z1 - self.z0))

            nlines = 500000   #Number from the original integration, not this one
            nx_lines = int(np.sqrt(nlines)*(self.x1-self.x0)/(self.y1 - self.y0))
            ny_lines = int(nlines/nx_lines)


            self.xsl = np.linspace(self.xs[0], self.xs[-1], self.nx)
            self.ysl = np.linspace(self.ys[0], self.ys[-1], self.ny)

            x, y = np.meshgrid(self.xsl, self.ysl)
            z = 0.0*np.ones((np.shape(x)))
            surface = pv.StructuredGrid(x, y, z)

            p.background_color = "darkslategrey"

            p.subplot(0, 0)
            vmax = np.percentile(np.abs(self.surface_array),99)
            p.add_mesh(surface, scalars = self.surface_array, show_edges=False,cmap = 'berlin', clim = [-vmax, vmax])

            camera_position = p.camera.position

            posx, posy, posz = camera_position
            xcentre =  (self.x1 + self.x0)/2
            ycentre =  (self.y1 + self.y0)/2

            #That doesn't work because angles. Bugger.
            print(posx, xcentre, posy, ycentre)
            xyabs = np.sqrt((posx-xcentre)**2 + (posy-ycentre)**2)
            #xyabs = xyabs/4
            zoom = 1.5

            new_position = [xcentre + xyabs*np.sin(plot_angle)/zoom, ycentre -xyabs*np.cos(plot_angle)/zoom, posz/zoom]

            p.camera.position = new_position

            p.camera.focal_point = ((self.x1 + self.x0)/2,(self.y1 + self.y0)/2,self.z0 + (self.z1 - self.z0)/4)

            p.remove_scalar_bar()

            p.show(screenshot='plots/%s_lines.png' % input_fname[:-3], window_size = (2000,2000))
            p.close()

    fltrace = Fltrace(input_fname = input_fname, nlines = nlines, show = show, id = id, do_flh = False, remove_files = remove_emission_files, swapaxes = swapaxes)
    if not justplot:
        fltrace.trace_lines()

    fltrace.plot_fieldlines()














