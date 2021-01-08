[![DOI](https://zenodo.org/badge/238084191.svg)](https://zenodo.org/badge/latestdoi/238084191)

# Learning dominant physical processes with data-driven balance models

Code and details for the [manuscript](https://arxiv.org/abs/2001.10019) on data-driven dominant balance analysis.  The repository should have everything needed to reproduce the results of the paper and get started exploring dominant balance in other systems.

The method essentially uses out of the box `scikit-learn` algorithms, so this repository doesn't have a separate library for import.  All the Jupyter notebooks should basically be self-contained, provided the dependencies (below) are met. See Appendix A of the paper for detailed information on data provenance, but there are notebooks for the following examples:

### 1. Burgers equation

Simple PDE model for shock formation.  This notebook is the most heavily annotated and is the easiest to understand and visualize, since the equation only has three terms.  Data was generated with second-order finite difference solver

### 2. $\mathrm{Re}=100$ cylinder wake

Vortex shedding past a circular cylinder at moderate Reynolds number. This notebook also includes a demonstration of uncertainty analysis for the dominant balance procedure.  Data was generated with a direct numerical simulation using the [Nek5000](https://nek5000.mcs.anl.gov/) spectral element solver.  The case files necessary to re-run the DNS are in the data folder.

### 3. Turbulent boundary layer

Main example from the paper.  In addition to the results in the paper, the notebook shows that the identified balance regions are consistent with two kinds of scaling collapse: the Blasius profile in the laminar inflow region, and wall coordinates in the viscous sublayer.  Data available from the [Johns Hopkins Turbulence Database](http://turbulence.pha.jhu.edu/Transition_bl.aspx). We only make use of the mean profile, which is available as an HDF5 file.

### 4. Supercontinuum generation

Optical pulse propagation with the Generalized Nonlinear Schrodinger Equation.  MATLAB spectral solver code originally obtained [here](www.scgbook.info), but modified scripts to run the solver, nondimensionalize the results, and extract the various terms from the solver are in the data folder.

### 5. Geostrophic balance

Surface currents in the Gulf of Mexico.  Data was downloaded from the [HYCOM](https://www.hycom.org/) group, but a MATLAB script to do the same is in the data folder.

### 6. Neuron model

Generalized Hodgkin-Huxley model of an intrinsically bursting neuron.  Again, MATLAB scripts implementing the full model are in the data folder.

### 7. Rotating detonation engine analogue

Modified Majda-Burgers model of 1D combustion in a rotating detonation engine

__Dependencies:__

* `scikit-learn`: all machine learning methods are implemented as off-the-shelf models
* `pymech`: python [package](https://github.com/jcanton/pymech) helpful for reading/writing Nek5000 data files (only needed for cylinder example)
* `pyclaw`: python interface to the [Clawpack](https://www.clawpack.org/) finite volume solver (only needed to generate data for the RDE model)
* Standard libraries:`seaborn`, `scipy`, `numpy`, `matplotlib`
