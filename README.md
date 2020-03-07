# PaCS-MD simulation script

* PaCS-MD-dissociation.py: python script for PaCS-MD dissociation simulations.
* PaCS-MD-docking.py: python script for switching between association and dissociation PaCS-MD simulations 


## Prerequisites
* GROMACS
* python 3.x with standard package and numpy.

```
This script is tested with GROMACS version 5.x.

### Known bug issues:

1. GROMACS 2019.x ++ on AMD gpu card: add -pme cpu in gpuid.
Example: gpuid="0 -pme cpu" and "gpu=1"

```
