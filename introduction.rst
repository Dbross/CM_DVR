:Author: David Bross <dbross@anl.gov> $
:Date: 2022-06-01 $
:Revision: 1 $
:Description: Multidimensional DVR program to generate and solve 1-d and 2-d potentials.

Introduction
************
This program was developed to calculate rovibrational levels by solving n-dimensional potentials. The program is capable of generating the necessary grid and, once the potential energy for those coordinates has been calculated, solving it for the given dimensionality. The potential is solved using the ‘universal’ Discrete Variable Representation (DVR) originally developed by Colbert and Miller [Colbert1992]. While the implementation can entertain both aperiodic and periodic potentials, the intrinsic target are soft and large-amplitude motions of non-rigid molecules. For the present purpose, non-rigid molecules are defined as those with at least one internal coordinate that undergoes a large amplitude motion [Bunker1998]. The potentials of these vibrational motions are poorly approximated with harmonic oscillators. The current implementation is in practice limited to solving 1 and 2 dimensional potentials due to the strongly scaling cost of these methods with increased dimensionality.


.. [Colbert1992] Colbert, D. T.; Miller, W. H. A Novel Discrete Variable Representation for Quantum-Mechanical Reactive Scattering via the S-Matrix Kohn Method. J. Chem. Phys. 1992, 96, 1982–1991. http://doi.org/10.1063/1.462100 
.. [Bunker1998] Bunker, P. R.; Jensen P. Molecular Symmetry and Spectroscopy, 2nd ed.; NRC Research Press: Ottawa 1998.

This code was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, Division of Chemical Sciences, Geosciences and Biosciences, under Contract No. DE-AC02-06CH11357, through the Computational Chemical Sciences Program (D.H.B.) and the Gas-Phase Chemical Physics Program (B.R.). 
