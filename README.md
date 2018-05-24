# GRID-core

GRID-core is a core-finding method using the contours of the local gravitational potential to identify core boundaries, 
as described in Gong & Ostriker (2011). This method can be applied to both observed 2D surface density and simulated 3D 
volume density. Detailed discussion and user's manual can be found here: http://www.astro.umd.edu/~hgong/GRID_core.htm 

An updated 2D version of the GRID-core algorithm in IDL is avilable here. This script is suitable for identifying 
clumps/cores in observed column density maps. The required input is a two-dimensional FITS file containing a map of the 
column density of atomic hydrogen. The code will generate two fits files containing information of identified cores, as 
well as a map illustrating the gravitational potential contours and the identified cores. A sample FITS map is provided 
here for code-testing purpose.

NOTE: This code uses publically-available IDL routines imdisp.pro and extremes_2d.pro, both included here for convenience.
