# EMSIPON
Mesoscopic Simulations of Polymers

Source code of a mesoscopic, mixed particle- and field-based Brownian dynamics methodology for the simulation of 
entangled polymer melts. 

Polymeric beads consist of several Kuhn segments, and their motion is dictated by the Helmholtz energy of the 
sample, which is a sum of the entropic elasticity of chain strands between beads, slip springs, and nonbonded 
interactions. 
 
The entanglement effect is introduced by the slip springs, which are springs connecting either nonsuccessive 
beads on the same chain or beads on different polymer chains. The terminal positions of slip springs are altered 
during the simulation through a kinetic Monte Carlo hopping scheme, with rate-controlled creation/destruction 
processes for the slip springs at chain ends. 

The rate constants are consistent with the free energy function employed and satisfy microscopic reversibility at 
equilibrium. The free energy of nonbonded interactions is derived from an appropriate equation of state, and it is 
computed as a functional of the local density by passing an orthogonal grid through the simulation box; accounting 
for it is necessary for reproducing the correct compressibility of the polymeric material. 

The mesoscopic simulation methodology has been implemented for the case of cis-1,4-polyisoprene, whose structure, 
dynamics, thermodynamics, and linear rheology in the melt state are quantitatively predicted and validated without 
a posteriori fitting the results to experimental measurements.

References

If you plan to publish an academic paper using this software, please consider citing the following 
publications:
1. Vogiatzis, G.G.; Megariotis, G.; Theodorou, D.N. Macromolecules 2017, 50, 3004. 
2. Megariotis G.; Vogiatzis, G.G.; Sgouros, A.P.; Theodorou D.N. Polymers 2019, 10, 1156.
