# fltrace
Basic cartesian magnetic field line tracer, for use with staggered grids.

This repository is set up to run using the example data in the folder '/Data/', which represents the initial condition for a solar jet model.

The code runs using the python wrapper 'fltrace_fast.py', which establishes the starts of each field line (in the function set_starts). It then generates arrays to allow variables to be read-in to fortran at run time, which are saved to the folder 'fl_tools'. 

# Steps to run:
Check the Makefile is pointing to the correct modules. This one uses gfortran and netcdf -- the locations of these are set at the top of the file. Shouldn't need to touch anything else below here.

Change the data directory paths in 'fltrace_fast.py' and 'src/fltrace.f90' if necessary. Currently set up to use '/Data/'. 

Run by calling the python script. The code will compile if anything has changed, but it's proably helpful to try doing so before running the python.

# Data information:
The magnetic field data are stored as netcdf files with the tags 'bx', 'by' and 'bz'. Take care with the order of the dimensions -- to read in to python usually the x and z dimensions need to be flipped relative to fortran. It's easier to flip in python than fortran so that's what the code does. 

The dimensions are the number of CELLS fully within the grid domain, not the number of points. One dimension of each magnetic field array is stored on the grid points due to the staggering -- hence it has a dimension one unit larger. 

