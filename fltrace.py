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
'''
id = 4        #Identifier for this run, so several can be done concurrently or in sequence

#Copy files to plot into input folder:
os.system('scp -r /extra/tmp/trcn27/mf3d/%03d/1000.nc ./input/1000.nc' % id)

#Define parameters for plotting:
#input_fname = 'Bout__2021.0621.061200.nc'    #File name within the directory ./input/'
input_fname = '1000.nc'    #File name within the directory ./input/'

nlines = 1000000                            #Approx number of field lines to trace
show = True                                #Plots various things, including the imported magnetic field on the lower boundary
justplot = False                          #If true, finds existing data and just plots it. If you want to just tweak the plots without running everything again.
closed_boundaries = True               #If there are field lines through the x and y boundaries, set to False. This does Chris' correction to the vector potentials, which is quite slow.
remove_emission_files = False               #Removes the emission files after plotting
swapaxes = True                             #Swaps x and z axes of the imported file

class Fltrace():
    def __init__(self, input_fname, id = 0, nlines = 10000, show = True, do_flh = False, closed_boundaries = False, swapaxes = False, remove_files = False):

        #Some parameters
        self.max_line_length = 100000
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

        if self.show: #Plot slice of each coordinate to check this looks right. If it doesn't, try swapping the import axes?
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

        if do_flh:
            print('Finding various vector potentials for helicity calculations...')
            flh = FLH(self, closed_boundaries = closed_boundaries)
            print('Seemingly complete.')
            self.flh_density = flh.flh_density
            self.winding_density = flh.winding_density
            self.twist_density = flh.twist_density
        print('Saving for use by Fortran')
        save_for_fortran(self)
        print('Saved to temporary file')

    def trace_lines(self):

        #Folder admin
        if not os.path.exists('./fl_data/'):
            os.mkdir('fl_data')

        #Find start points. Assumed to be on a regularly-spaced grid
        self.set_starts()

        #Create runtime variables for fortran and save out to tmp folder
        self.setup_tracer()
        #Do the tracing
        self.trace_lines_fortran()

        #Remove temporary files
        os.system('rm ./fl_data/flparameters%05d.txt' % self.id)
        os.system('rm ./fl_data/starts%05d.txt' % self.id)
        os.system('rm ./tmp/%05d.nc' % (self.id))



    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.

        #Plot up from the surface, based on some threshold of how strong the magnetic field is...
        nlines = self.nlines

        nx_lines = int(np.sqrt(nlines)*(self.x1-self.x0)/(self.y1 - self.y0))
        ny_lines = int(nlines/nx_lines)

        alpha = 2.0

        self.starts = []

        boundary_nolines = 0.9  #Fraction of the domain to actually bother plotting from
        dx_all = self.xs[-1] - self.xs[0]
        dy_all = self.ys[-1] - self.ys[0]
        xstart = self.xs[0] + 0.5*(1.0 - boundary_nolines)*dx_all
        xend = self.xs[-1] - 0.5*(1.0 - boundary_nolines)*dx_all
        ystart = self.ys[0] + 0.5*(1.0 - boundary_nolines)*dy_all
        yend = self.ys[-1] - 0.5*(1.0 - boundary_nolines)*dy_all
        x_allstarts = np.linspace(xstart, xend, nx_lines)
        y_allstarts = np.linspace(ystart, yend, ny_lines)
        for i in range(nx_lines):
             for j in range(ny_lines):
                 #This only works if symmetric!
                 xstart = x_allstarts[i]
                 ystart = y_allstarts[j]

                 self.starts.append([xstart, ystart, self.zs[0]])

        self.starts = np.array(self.starts)
        fig = plt.figure(figsize = (5,5))
        plt.scatter(self.starts[:,0], self.starts[:,1], s= 0.1, c= 'black')
        plt.close()
        print('Tracing', len(self.starts), 'lines')

        self.nstarts = len(self.starts)
        self.starts = self.starts.reshape(self.nstarts*3)

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.1 #Tracing 'timestep' as a proportion of the self size
        print_flag = 1  #Print some things as the tracing happens
        save_all = 0

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
        if os.uname()[1] == 'brillouin.dur.ac.uk' or os.uname()[1] == 'modigliani.dur.ac.uk':
            os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace %d' % self.id)
        elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
            os.system('mpiexec -np 1 ./bin/fltrace %d' % self.id)
        else:
            os.system('mpirun -np 1 ./bin/fltrace %d' % self.id)

    def plot_emissions(self, allscales = [], remove_files = True):

        try:
            data = netcdf_file('./fl_data/emiss%05d.nc' % self.id, 'r', mmap=False)
            print('Emissions found with id', self.id)

        except:
            print('File not found -- tracing has probably failed. Sorry.')

        self.xsum = np.swapaxes(data.variables['emiss_xsum'][:],0,1)
        self.ysum = np.swapaxes(data.variables['emiss_ysum'][:],0,1)
        self.zsum = np.swapaxes(data.variables['emiss_zsum'][:],0,1)

        self.flh_array = data.variables['flh_array'][:]
        self.winding_array = data.variables['winding_array'][:]
        self.twist_array = data.variables['twist_array'][:]
        self.surface_array = data.variables['surface_array'][:]
        self.current_array = data.variables['current_array'][:]

        data.close()

        #Using the start points already established, set this array up
        nlines = self.nlines
        nx_lines = int(np.sqrt(nlines)*(self.x1-self.x0)/(self.y1 - self.y0))
        ny_lines = int(nlines/nx_lines)

        flh_mesh = np.zeros((nx_lines, ny_lines))
        winding_mesh = np.zeros((nx_lines, ny_lines))
        twist_mesh = np.zeros((nx_lines, ny_lines))
        surface_mesh = np.zeros((nx_lines, ny_lines))
        current_mesh = np.zeros((nx_lines, ny_lines))

        count = 0
        for i in range(nx_lines):
             for j in range(ny_lines):
                flh_mesh[i,j] = self.flh_array[count]
                winding_mesh[i,j] = self.winding_array[count]
                twist_mesh[i,j] = self.twist_array[count]
                surface_mesh[i,j] = self.surface_array[count]
                current_mesh[i,j] = self.current_array[count]
                count += 1

        xplots = np.linspace(self.xs[0], self.xs[-1], nx_lines)
        yplots = np.linspace(self.ys[0], self.ys[-1], ny_lines)

        fig, axs = plt.subplots(2,4, figsize = (10, 7))

        ax = axs[0,0]
        toplot = self.bz[:,:,0].T
        ax.set_title("Magnetogram")
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(self.xs, self.ys, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_aspect('equal')

        ax = axs[0,1]
        toplot = flh_mesh.T
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Field-line helicity")
        ax.set_aspect('equal')

        ax = axs[1,1]
        toplot = flh_mesh.T*np.abs(surface_mesh.T)
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Weighted field-line helicity")
        ax.set_aspect('equal')

        ax = axs[1,0]
        toplot = winding_mesh.T
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Winding")
        ax.set_aspect('equal')

        ax = axs[0,2]
        toplot = twist_mesh.T
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Twist")
        ax.set_aspect('equal')

        ax = axs[1,2]
        toplot = twist_mesh.T*np.abs(surface_mesh.T)
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Weighted twist")
        ax.set_aspect('equal')

        ax = axs[0,3]
        toplot = current_mesh.T
        vmax = np.percentile(np.abs(toplot), 99)
        vmin = -vmax
        ax.pcolormesh(xplots, yplots, toplot,cmap ='seismic', vmax = vmax, vmin = vmin)
        ax.set_title("Total Current along line")
        ax.set_aspect('equal')

        ax = axs[1,3]
        ax.set_axis_off()

        plt.tight_layout()
        plt.savefig('./plots/%s_flh.png' % self.input_fname[:-3])
        plt.show()

        xsplot = np.linspace(self.x0,self.x1,self.zsum.shape[0]+1)
        ysplot = np.linspace(self.y0,self.y1,self.zsum.shape[1]+1)
        zsplot = np.linspace(self.z0,self.z1,self.xsum.shape[1]+1)

        #fig, axs = plt.subplots(4,1, figsize = (15,7.5))
        fig = plt.figure(figsize = (10,7.5))

        toplots = [np.sqrt(self.zsum).T, np.sqrt(self.ysum), np.sqrt(self.xsum)]
        ranges = [[ysplot, xsplot], [xsplot, zsplot],[ysplot, zsplot]]
        titles = ['Top', 'x face', 'y face']

        gs0 = gridspec.GridSpec(4,2,figure=  fig)


        for i in range(3):
            if i == 0:
                ax = fig.add_subplot(gs0[0:2,0])
            else:
                ax = fig.add_subplot(gs0[2:4,i-1])
            toplot = toplots[i]

            #Top down
            #vmax = np.percentile(toplot, 99.75)
            pc = (np.sum([toplot > 1e-5])/np.sum([toplot >= 0.0]))
            pc_total = 100 - 0.1*pc
            vmax = np.percentile(toplot, pc_total)
            vmax = max(0.1, vmax)

            if len(allscales) > 0:
                vmax = allscales[self.snap, i]
            #vmax = 1.0
            im = ax.pcolormesh(ranges[i][0], ranges[i][1],toplot.T, vmin = 0, vmax = vmax, cmap = 'inferno')
            ax.set_aspect('equal')
            ax.set_title(titles[i])


        #Then magnetogram
        ax = fig.add_subplot(gs0[0:2,1])
        toplot = self.bz[:,:,0]
        ax.pcolormesh(self.ys, self.xs, toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
        ax.set_aspect('equal')
        ax.set_title('Lower boundary magnetogram')

        plt.tight_layout()
        plt.savefig('./plots/%s_whitelight.png' % self.input_fname[:-3])
        if self.show:
            plt.show()
        plt.close()

        if self.remove_files:
            os.system('rm -r ./fl_data/emiss%05d.nc' % self.id)

fltrace = Fltrace(input_fname = input_fname, nlines = nlines, show = show, id = id, do_flh = not justplot, closed_boundaries = closed_boundaries, remove_files = remove_emission_files, swapaxes = swapaxes)
if not justplot:
    fltrace.trace_lines()
fltrace.plot_emissions(remove_files = False)















