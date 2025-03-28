Field line tracer and tools to calculate field line helicity, winding and twist. Currently set up for use with solar active regions.

To run: 

1. Put input files in a folder called ./inputs/. Or change the input folder in fltrace.py
2. Set up parameters at the top of fltrace.py. Options to just plot or number of field lines to calculate etc. are given.
3. Dependencies should all be obvious *except* for if you'd like to export the file as an interactive html. To do that you'll need

```
pip install pyvista trame
pip install 'pyvista[jupyter]'
```

4. IMPORTANT - change the locations of gfortran and netcdf at the top of the Makefile. Currently set up for a linux machine in the maths department.
5. Run just by running the wrapper fltrace.py

Plots should pop up or be saved to ./plots/

Dependencies:

numpy
scipy
matplotlib
pyvista
streamtracer

I think that's all
