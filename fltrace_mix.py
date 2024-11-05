import numpy as np
import matplotlib.pyplot as plt
import random
import os
from scipy.io import netcdf_file
import fltools

import pyvista as pv
pv.start_xvfb()

#Compile fortran:
#os.system('f2py3 -c -m fltools fltools.f90')

class trace_fieldlines():
    def __init__(self):
        #Establish grid parameters (can be read in from elsewhere of course)
        self.run = 0
        self.snap = 0
        self.print_flag = 1
        self.save_number = self.snap
        self.data_root = './Data/'


        data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)

        bx_import = np.swapaxes(data.variables['bx'][:],0,2)
        by_import = np.swapaxes(data.variables['by'][:],0,2)
        bz_import = np.swapaxes(data.variables['bz'][:],0,2)

        self.nx = np.shape(bz_import)[0]
        self.ny = np.shape(bz_import)[1]
        self.nz = np.shape(bz_import)[2] - 1

        self.bx = np.zeros((self.nx+1,self.ny+2,self.nz+2))
        self.by = np.zeros((self.nx+2,self.ny+1,self.nz+2))
        self.bz = np.zeros((self.nx+2,self.ny+2,self.nz+1))

        self.bx[:,1:-1,1:-1] = bx_import
        self.by[1:-1,:,1:-1] = by_import
        self.bz[1:-1,1:-1,:] = bz_import
        data.close()

        self.x0 = -12.; self.x1 = 12.
        self.y0 = -12.; self.y1 = 12.
        self.z0 = -24.0/self.nz; self.z1 = 24

        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        self.dz = self.zs[1] - self.zs[0]

        #Establish start points for the field line plotting
        self.max_line_length = 10000
        self.ds = 0.1*min(self.dx, self.dy, self.dz) #Tracing 'timestep' as a proportion of the grid size
        self.weakness_limit = 1e-3   #Minimum field strength to stop plotting
        self.line_plot_length = 50  #To save time while plotting, reduce the length of the plotted lines

        #Folder admin
        if not os.path.exists('./fl_data/'):
            os.mkdir('fl_data')
        os.system('rm ./fl_data/flines.nc')
        #Find start points
        self.set_starts()
        #Create runtime variables for fortran
        #self.setup_tracer()
        #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
        self.trace_lines_python()
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
        p.add_mesh(surface, scalars= self.bz[1:-1,1:-1,0], show_edges=True,cmap = 'plasma')

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

    def trace_lines_python(self):

        self.lines = np.zeros((self.nstarts*2, self.max_line_length, 3))
        for li, start in enumerate(self.starts[:]):
            self.lines[2*li] = fltools.integrate_line(self.bx,self.by,self.bz,start,self.max_line_length,self.nx,self.ny,self.nz,self.x0,self.x1,self.y0,self.y1,self.z0,self.z1,self.ds,self.weakness_limit,-1)
            self.lines[2*li+1] = fltools.integrate_line(self.bx,self.by,self.bz,start,self.max_line_length,self.nx,self.ny,self.nz,self.x0,self.x1,self.y0,self.y1,self.z0,self.z1,self.ds,self.weakness_limit,1)

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
        self.starts = np.array(self.starts)

    def interpolate_field(self,pt):
        #Outputs the magnetic field vector at point pt, using the individual magnetic field vectors (no unecessary averaging)
        #Find position of point on the grid
        b1 = np.zeros((3))
        xp = (self.nx)*(pt[0] - self.x0)/(self.x1 - self.x0)
        yp = (self.ny)*(pt[1] - self.y0)/(self.y1 - self.y0)
        zp = (self.nz)*(pt[2] - self.z0)/(self.z1 - self.z0)

        #Interpolate bx
        xi = int(xp); yi = int(yp + 0.5); zi = int(zp + 0.5)
        xf = xp - xi; yf = yp+0.5 - yi; zf = zp + 0.5 - zi

        b1[0] = b1[0] + self.bx[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.bx[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.bx[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.bx[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[0] = b1[0] + self.bx[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.bx[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.bx[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.bx[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        #Interpolate by
        xi = int(xp+0.5); yi = int(yp); zi = int(zp + 0.5)
        xf = xp + 0.5 - xi; yf = yp - yi; zf = zp + 0.5 - zi
        b1[1] = b1[1] + self.by[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.by[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.by[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.by[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[1] = b1[1] + self.by[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.by[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.by[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.by[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        #Interpolate bz
        xi = int(xp+0.5); yi = int(yp + 0.5); zi = int(zp)
        xf = xp + 0.5 - xi; yf = yp + 0.5 - yi; zf = zp - zi
        b1[2] = b1[2] + self.bz[xi,yi,zi]*(1-xf)*(1-yf)*(1-zf) + self.bz[xi,yi,zi+1]*(1-xf)*(1-yf)*(zf) + self.bz[xi,yi+1,zi]*(1-xf)*(yf)*(1-zf) + self.bz[xi,yi+1,zi+1]*(1-xf)*(yf)*(zf)
        b1[2] = b1[2] + self.bz[xi+1,yi,zi]*(xf)*(1-yf)*(1-zf) + self.bz[xi+1,yi,zi+1]*(xf)*(1-yf)*(zf) + self.bz[xi+1,yi+1,zi]*(xf)*(yf)*(1-zf) + self.bz[xi+1,yi+1,zi+1]*(xf)*(yf)*(zf)

        return (b1/np.sqrt(np.sum(b1**2)), np.sqrt(np.sum(b1**2)))

    def trace_line(self, start, updown = 0):
        #Traces an individual field line
        pts = 1e6*np.ones((self.max_line_length, 3))
        count = 0
        go = True
        pts[0] = start
        pt = start
        if not self.inbounds(pt):
            go = False
        else:
            grad, _ = self.interpolate_field(pt)
            mag = np.sqrt(np.sum(grad**2))
            grad = grad/mag

        while go:
            if not self.inbounds(pt):
                go = False
            elif mag < self.weakness_limit:
                go = False
            elif count > self.max_line_length-1:
                go = False
            else:
                #grad, _ = self.interpolate_field(pt)
                grad = fltools.interpolate_field(self.bx, self.by, self.bz, pt, self.nx,self.ny,self.nz,self.x0,self.x1,self.y0,self.y1,self.z0,self.z1)
                mag = np.sqrt(np.sum(grad**2))
                grad = grad/mag
                pt = pt + updown*self.ds*grad
                pts[count] = pt
                count = count + 1
        return pts

    def inbounds(self, pt):
        #returns true if the point is still within the domain. Otherwise doesn't
        if pt[0] < self.x0 or pt[0] > self.x1 or pt[1] < self.y0 or pt[1] > self.y1 or pt[2] < self.z0 or pt[2] > self.z1:
            return False
        else:
            return True

trace_fieldlines()


