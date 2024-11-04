import numpy as np
import matplotlib.pyplot as plt
import random
import os
from scipy.io import netcdf_file

if os.uname()[1] == 'brillouin.dur.ac.uk':
    import pyvista as pv
    pv.start_xvfb()


class trace_fieldlines():
    def __init__(self, grid, bx,by,bz,save = -1,plot_vista = True,plot_notvista=False):
        print('Tracing field lines...')
        self.xs = np.linspace(grid.x0,grid.x1,grid.nx+1)
        self.ys = np.linspace(grid.y0,grid.y1,grid.ny+1)
        self.zs = np.linspace(grid.z0,grid.z1,grid.nz+1)

        self.nx = grid.nx
        self.ny = grid.ny
        self.nz = grid.nz

        self.x0 = grid.x0; self.x1 = grid.x1
        self.y0 = grid.y0; self.y1 = grid.y1
        self.z0 = grid.z0; self.z1 = grid.z1

        self.xc = np.zeros(self.nx + 2)
        self.yc = np.zeros(self.ny + 2)
        self.zc = np.zeros(self.nz + 2)

        self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
        self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
        self.zc[1:-1] = 0.5*(self.zs[1:] + self.zs[:-1])

        self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
        self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])
        self.zc[0] = self.zc[0] - (self.zc[2] - self.zc[2])

        self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
        self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])
        self.zc[-1] = self.zc[-2] + (self.zc[-2] - self.zc[-3])

        self.dx = np.sum(self.xs[1:] - self.xs[:-1])/len(self.xs[1:])
        self.dy = np.sum(self.ys[1:] - self.ys[:-1])/len(self.ys[1:])
        self.dz = np.sum(self.zs[1:] - self.zs[:-1])/len(self.zs[1:])

        self.bx = bx; self.by = by; self.bz = bz
        self.ds = 0.1*min(self.dx,self.dy,self.dz)

        self.save = save
        #Find start positions of each field line. For now a regular grid on the lower boundary
        if False: #Cartesian tracing grid
            nxs = 10; nys = 12
            self.starts = []

            xis = np.linspace(self.x0+1e-6,self.x1-1e-6,nxs)
            yjs = np.linspace(self.y0+1e-6,self.y1-1e-6,nys)

            for i in range(nxs):
                for j in range(nys):
                    self.starts.append([xis[i],yjs[j],1e-6])
        else:   #Polar plotting grid
            self.starts = []
            #Trace from the top
            nrs = 30; nthetas = 2
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

        #self.starts = [[0,10.0,0]]

        if False:
            self.lines = []
            for start in self.starts:
                doplot = True
                line = self.trace_line(start)
                if len(line) < 1:
                    continue

                if line[-1][2] > 1.0:
                    doplot = False
                if doplot:
                    self.lines.append(line)

        else:
            self.setup_tracer()
            self.trace_lines_fortran()

        #Do use pyvista
        x, y = np.meshgrid(self.xs, self.ys)
        z = 0*x*y
        surface = pv.StructuredGrid(x, y, z)
        p = pv.Plotter(off_screen=True)
        p.background_color = "black"
        p.add_mesh(surface, scalars= bz[1:-1,1:-1,0], show_edges=True,cmap = 'plasma')

        for li, line in enumerate(self.lines):

            line = np.array(line)
            line_length = len(line[line[:,2]<1e6])

            line = np.array(line[:line_length]).tolist()

            doplot = True
            if line_length == 0:
                doplot = False
            elif line[0][2] < 1.0 and line[-1][2] > 20.0:
                doplot = False

            if doplot:

                p.add_mesh(pv.Spline(line, len(line)),color='white')

        p.camera.position = (20.0,40,20.0)
        p.camera.focal_point = (0,0,4)
        p.show(screenshot='plots/b%04d.png' % save, window_size = (1000,1000))


    def setup_tracer(self):

        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.05 #Tracing 'timestep' as a proportion of the grid size
        weakness_limit = 1e-3   #Minimum field strength to stop plotting
        print_flag = 1  #Print some things as the tracing happens

        nstarts = len(self.starts)
        starts = np.array(self.starts).reshape(nstarts*3)

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
        variables[10] = self.save
        variables[11] = nstarts
        variables[12] = print_flag
        variables[13] = max_line_length
        variables[14] = ds
        variables[15] = weakness_limit

        np.savetxt('flparameters.txt', variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('starts.txt', starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):
        os.system('make')
        os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace')
        data = netcdf_file('flines.nc', 'r', mmap=False)

        try:
            data = netcdf_file('flines.nc', 'r', mmap=False)
            print('Field lines found')

        except:
            print('File not found')

        self.lines = np.swapaxes(data.variables['lines'][:],0,2)

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

        return b1/np.sqrt(np.sum(b1**2))

    def trace_line(self, start):
        #Traces an individual field line
        go = True
        pts = [[start[0],start[1],start[2]]]
        pt = start
        if not self.inbounds(pt):
            return []
        grad = self.interpolate_field(pts[-1])
        test = pt + self.ds*grad
        if self.inbounds(test):
            updown = 1
        else:
            updown = -1

        while go:
            if not self.inbounds(pt):
                go = False
            else:
                grad = self.interpolate_field(pts[-1])
                pt = pt + updown*self.ds*grad
                pts.append([pt[0],pt[1],pt[2]])

        return pts

    def inbounds(self, pt):
        #returns true if the point is still within the domain. Otherwise doesn't
        if pt[0] < self.x0 or pt[0] > self.x1 or pt[1] < self.y0 or pt[1] > self.y1 or pt[2] < self.z0 or pt[2] > self.z1:
            return False
        else:
            return True



