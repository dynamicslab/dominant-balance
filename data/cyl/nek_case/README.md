# Flow past a cylinder in 2D

Modified from `ext_cyl` example case. 

## gmsh2nek
The geometry is meshed with gmsh and has to be converted to the .re2 format. The gmsh2nek utility works, but the settings have to be as follows:

* 2D mesh (so not like the pseudo-2D for OpenFOAM)
* Regions should be labeled with "Physical Line" in .geo file
* Set order 2 in options
* Export mesh to .msh file with version 2
* Run gmsh2nek to convert to .re2 format (note this needs a .par parameter file now) and genmap to parallelize.
* Set boundary conditions in usrdat2 subroutine in the .usr file. 