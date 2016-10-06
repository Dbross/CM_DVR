:Author: David Bross <dbross@anl.gov>
:Date: 2016-10-06 $
:Revision: 1 $
:Description: Multidimensional DVR program to generate and solve n-d potentials.

Introduction
************
This program was developed to solve n-dimensional potentials grids. Its was developed to calculate frequencies of soft and non-rigid molecular motions. It is capable of generating the necessary grid and once the potential energy for these coordinates has been calculated, solving that grid. The potential is solved with an implementation of the Kohn variational method using  Discrete Variable Representation originally developed by Colbert and Miller [Colbert1992]_.
Non-rigid molecules are defined the same as in [Bunker]_ meaning those with at least one internal coordinates that undergoes a large amplitude motion. A large amplitude motion is a change in the internal coordinate that is of the same magnitude as the internal coordinate itself (e.g. a change of equal to the bond length or angle of the internal coordinate). The potentials of these vibrational motions are poorly approximated with harmonic fits of the potential.


.. [Colbert1992] Colbert, D. T., & Miller, W. H. (1992). A Novel Discrete Variable Representation for Quantum-Mechanical Reactive Scattering via the S-Matrix Kohn Method. The Journal of Chemical Physics, 96(3), 1982–1991. http://doi.org/10.1063/1.462100
.. [Bunker] Bunker's book...
