# FQparam
script useful for the parametrization of FQ polarizable force field

Goal: optimize the parameters to minimize the difference between interaction energies and multipoles of a QM molecules in the presence of a fixed charge and of the same molecules with FQ charges. The fized charges used for the parametrization are chosen to be the centers of PCM tesserea but it can be modified.

Structure:  
Input:  
* coordinates
* atom types
* Definition of functional to be minimized
* path for a database with charge positions, intensity, QM/MM interaction energy, dipole and quadrupole of QM density


Scripts:  
* get position of charges from PCM cavity
* prepare openmolcas input file following a template calculation
* parse Openmolcas log and retrieve: charge position and intesity, Interaction energy, dipole and quadrupole
* Optimizer

Optimizer:  
read input file, for each charge position compute and store:
* FQ charges
* FQ interaction energy
* FQ dipole and quadrupole
* derivative of FQ charges with respect to the parameters
* compose the derivative of the selected functional
* use optimization procedure to find minimum
* compute statistics on the optimization

NB: use reasonable measure unit for energy interaction and for multipoles
