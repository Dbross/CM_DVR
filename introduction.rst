:Author: David Bross <dbross@anl.gov> $
:Date: 2019-04-02 $
:Revision: 1 $
:Description: Multidimensional DVR program to generate and solve 1-d and 2-d potentials.

Introduction
************
This program was developed to solve n-dimensional potentials grids. Its was developed to calculate frequencies of soft and non-rigid molecular motions. It is capable of generating the necessary grid and once the potential energy for these coordinates has been calculated, solving that grid. The potential is solved with an implementation of the Kohn variational method using  Discrete Variable Representation originally developed by Colbert and Miller [Colbert1992]_.
Non-rigid molecules are defined the same as in [Bunker]_ meaning those with at least one internal coordinates that undergoes a large amplitude motion. A large amplitude motion is a change in the internal coordinate that is of the same magnitude as the internal coordinate itself (e.g. a change of equal to the bond length or angle of the internal coordinate). The potentials of these vibrational motions are poorly approximated with harmonic fits of the potential.
The current implementation is limited to solving 1 and 2 dimensional potentials due to the strongly scaling cost of these methods with increased dimensionality.


.. [Colbert1992] Colbert, D. T.;  Miller, W. H. A Novel Discrete Variable Representation for Quantum-Mechanical Reactive Scattering via the S-Matrix Kohn Method. J. Chem. Phys. **1992**, 96, 1982–1991. http://doi.org/10.1063/1.462100
.. [Bunker] Bunker, P. R.; Jensen P. *Molecular Symmetry and Spectroscopy*, 2nd ed.; NRC Research Press: Ottawa 1998.
