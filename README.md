# PaCS-MD simulation script

* PaCS-MD-dissociation.py: python script for PaCS-MD dissociation simulations.
* PaCS-MDparallel.py: python script for PaCS-MD dissociation simulations performing parallelly all the replicas within a cycle parallelly. In addition, the analysis loop is also parallelized.
* dPaCS-MDparallelGROMACS.py: special version of the above for FUGAKU supercomputer using GROMACS
* adPaCSMD.py: python script for switching between association and dissociation PaCS-MD simulations 
* rmsdPaCSamber.py: python script for rmsd PaCS-MD using AMBER.


## Prerequisites
* GROMACS or AMBER 
* python 3.x with standard package and numpy.

```
This script is tested with GROMACS version 5.x.

### Known bug issues:

1. GROMACS 2019.x ++ on AMD gpu card: add -pme cpu in gpuid.
Example: gpuid="0 -pme cpu" and "gpu=1"

```
