/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @version 1.0 (Created on April 30, 2012)
 *  @brief The main C++ source file driving the execution of the code.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Georgios Vogiatzis and Grigorios Megariotis. 
 *  All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */


// Keep in mind the following web-page concerning Doxygen documentation system:
// https://modelingguru.nasa.gov/docs/DOC-1811

#include "b3D_integrator.h"
#include "Auxiliary.h"
#include "netmin.h"

using namespace NetworkNS;

/** @param[in] argc The count of arguments provided through the command line.
 *  @param[in] argv A 2D-array of characters containing the arguments provided through the 
 *                  command line.
 */
int main(int argc, char** argv) {

    
   unsigned int bd3D_nsteps = 0, bd3D_write_every = 0;
   
   cout.precision(12);
   
   cout << "#: Use of the code: "<< endl;
   cout << "#: ./netmin.x [data_file] [slip-spring creation rate] [nsteps] [nevery_steps]\n" << endl;
   
   std::string filename;
   if (argc > 1)
      filename = std::string(argv[1]);
   else{
      cout << "#: Data file has not been specified from command line.\n";
      filename = "emsipon.data";
      cout << "#: --- '" << filename << "' will be used as data file.\n";
   }
   
   double creation_rate;
   if (argc > 2)
      creation_rate = StringToNumber<double>(argv[2]);
   else
      return(1);
   
   if (argc > 3)
      bd3D_nsteps = atol(argv[3]);
   else
      bd3D_nsteps = 2000000000;
   
   if (argc > 4)
      bd3D_write_every = atol(argv[4]);
   else
      bd3D_write_every = 1000;
      

#ifdef FENE_SLS  
   cout << "#: (warning): FENE potential will be used for slip-springs.\n";
#endif
   
   // Define a new pointer to the class:
   NetwMin netw_min(filename);
   cb3D_integrator b3D_integrator(&netw_min, 500.0, creation_rate);
   
   cout << "#: Integrator has been initialized. Simulation will start for "<<bd3D_nsteps<< " steps.\n";
   cout << "#: The rate for the slip-spring creation is: " << creation_rate << endl;
   
   // Simulation starts here:
   double timestep = 1.0; // integration timestep in ps.
   cout << "#: The timestep of the integrator will be: " << timestep << " ps.\n";
   cout << "#: Report will be written every: " << bd3D_write_every << " steps.\n";
   
   
   b3D_integrator.integrate(bd3D_nsteps, timestep, bd3D_write_every);
     
   // Close trajectory file:
   netw_min.my_traj_file->~dump();
   
   
   return (EXIT_SUCCESS);
}

/** @mainpage Introduction
 *
 *  @section model Model 
 *  We have developed a Field Theory - inspired Monte Carlo approach for equilibrating long-chain
 *  melts based on a coarse-grained Hamiltonian that includes bonded, excluded volume, and cohesive
 *  interactions \cite Macromolecules_46_4670. Inputs to this method are the Kuhn length @f$ b @f$,
 *  the equilibrium melt density @f$ \rho_0 @f$ and the isothermal compressibility of the melt 
 *  @f$ \kappa_{\rm T} @f$ under the conditions of interest. Chains, which are represented as freely
 *  jointed sequences of Kuhn segments, adopt unperturbed conformations at equilibrium. This 
 *  scheme can be coarse-grained further by lumping sequences of @f$ n_{\rm Kuhns/bead} @f$  
 *  Kuhn segments into one bead.  
 *  The representation described above relies on soft intermolecular interactions between the beads 
 *  of different chains. These molecules can therefore cross each other, and do not exhibit any 
 *  signature of entanglement. Literature attempts to describe entanglements have largely relied on 
 *  so-called tube or network models. Our description of entanglements is based on the concept of 
 *  slip-springs @cite PhysRevLett_109_148302 .
 *  These slip-springs restrict the lateral motion of the beads with respect to the
 *  axis of the backbone, and favor a reptating motion along the backbone. 
 *  
 *  @section interactions Interactions
 *  At this coarse-grained level of representation, the polymeric system is assumed to be governed by a 
 *  Helmholtz free energy function @f$ A @f$ depending on its spatial extent, on the positions of 
 *  coarse-grained segments, on the connectivity, and on temperature. @f$ A @f$ is written as a sum of five
 *  terms reflecting contributions from the conformational entropy of strands between beads, slip-springs and
 *  nonbonded (polymer, polymer-particle and particle-particle) interactions
 *  @f[ A = A_{\rm b} + A_{\rm ss} + A_{\rm nb} + A_{\rm ps} + A_{\rm pp} @f]
 
 *  @subsection bonded_interactions Intramolecular Interactions
 *  
 *  The intramolecular interactions acting along the chain's backbone, between two beads @f$ i @f$ 
 *  and @f$ j @f$ located at positions @f$ \mathbf{R}_{i} @f$ and @f$ \mathbf{R}_j @f$, 
 *  respectively, are given by 
 *  @f[ A_{\rm b} \left(r^2_{ij} \right) = \epsilon_{\rm b} T 
 *  \frac{r^2_{ij}}{\sigma^2_{\rm b}}
 *  = \frac{3}{2}k_{\rm B}T\frac{\mathbf{R}_{ij} \cdot \mathbf{R}_{ij}}
 *  {n_{\rm Kuhns/bead} b^2} @f]
 *  with the parameters @f$ \epsilon_{\rm b} @f$ and @f$ \sigma_{\rm b} @f$ read from the 
 *  data file. In our approach @f$ \epsilon_{\rm b} = 0.012471 @f$  kJ/mol/K and 
 *  @f$ \sigma_{\rm b}^2 = 810 \; \text{\AA}^2 @f$ for polyisoprene melt.
 *  
 *  @subsection ss_interactions Slip-spring Interactions
 *  The slip-spring forcefield can be either modeled by an harmonic potential, similar to the one 
 *  used for the chain backbone, or by a Finitely Extensible Non-linear Elastic (FENE) potential of 
 *  the form 
 *  @f[ \mathcal{V}_{\rm ss} \left(r^2_{ij} \right) = - \epsilon_{\rm ss} T \sigma^2_{\rm ss} 
 *  \ln{\left[1- \frac{r^2_{ij}}{ \sigma^2_{\rm ss}} \right]} @f] 
 *  with the parameters @f$ \epsilon_{\rm ss}@f$ measured in kJ/mol/K/@f$\text{\AA}^2@f$ and 
 *  @f$ \sigma_{\rm ss} @f$ in @f$ \text{\AA}^2 @f$. Chappa et al. @cite PhysRevLett_109_148302, 
 *  give the following values to the above coefficients. The distance term can be connected to 
 *  our harmonic springs as  
 *  @f[\sigma^2_{\rm ss} = 5.76 \left( n_{\rm Kuhns/bead} b^2 \right) = 4665.6 \; \text{\AA}^2 @f]
 *  which is commensurate with the entanglement tube diameter of the polyisoprene.
 *  The energetic contribution, following the same authors, can be given as 
 *  @f[ \epsilon_{\rm ss} T \sigma^2_{\rm ss}  = \frac{k_{\rm ss}r^2_{\rm ss}}{2} 
 *  = 17.28 k_{\rm B} T @f]
 *  @f[\epsilon_{\rm ss} = \frac{17.28 k_{\rm B}}{\sigma^2_{\rm ss}} = 3.0794 \times 10^{-5} \;
 *  \frac{\rm kJ}{{\rm mol\:K}\:\text{\AA}^2} @f]  
 *  Alternatively, the strength of the slip-springs can be considered comparable to the strength 
 *  of the harmonic bonds along the chain backbone: 
 *  @f[ \epsilon_{\rm ss} T \sigma^2_{\rm ss}  = \frac{3}{2} k_{\rm B} T @f]
 *
 *  @subsection non_bonded_interactions Polymer Non-bonded Interactions
 *  To deal with nonbonded (excluded volume and van der Waals attractive) interactions in the 
 *  network representation, we introduce a network free energy: 
 *  @f[ A_{\rm nb} = \int d^3 \mathbf{r} f\left[\rho\left(\mathbf{r}\right) \right] @f]
 * 
 *  In the above equation, @f$\rho\left(\mathbf{r}\right) @f$ is the local density (number of
 *  Kuhn segments per unit volume) at position @f$ \mathbf{r} @f$ and 
 *  @f$ f\left(\rho \right) @f$ is a free energy density (free energy per unit volume). 
 *  Expressions for @f$f\left(\rho \right) @f$ may be extracted from an equation of state.
 *  Here the plan is to invoke a simple expression for @f$f\left(\rho \right) @f$, in the 
 *  form of a Taylor expansion, 
 *  @f[ f\left(\rho \right) = C\rho + B\rho^2 @f]
 *  with @f$ C @f$, @f$ B @f$ fitted such that the volumetric properties (pressure and 
 *  compressibility at mean density of interest) are reproduced.   
 *  Local density will be resolved only at the level of entire cells, defined by passing an 
 *  orthogonal grid through the entire system. The free energy of the system is approximated 
 *  by@f[ A_{\rm nb} = \sum_{\rm cells} V_{\rm cell}^{\rm acc} 
 *  f\left( \rho_{\rm cell}\right) @f]
 *  where @f$ V_{\rm cell}^{\rm acc} @f$, the accessible of a cell, is the volume of the 
 *  rectangular parallelepiped defining the cell minus the volume of any parts of 
 *  nanoparticles that may find themselves in the cell. 
 *  
 *  The cell density @f$\rho_{\rm cell} @f$ must be defined based on the nodal
 *  points in and around the cell, each nodal point contributing a mass equal to the node's 
 *  mass. Each nodal point @f$ j @f$ has mass @f$ n_j @f$ (in Kuhn segments) and a 
 *  characteristic size @f$ R_j @f$. We will discuss below how these quantities depend on the 
 *  node's molecular characteristics. We denote the position vector of node @f$ j @f$ by 
 *  @f$ \mathbf{r}_j = \left( x_j, y_j, z_j \right) @f$.
 *  The cell dimensions along the @f$ x @f$, @f$ y @f$, @f$ z @f$ directions will be denoted 
 *  as @f$ L_x @f$, @f$ L_y @f$, @f$ L_z @f$, respectively.
 *
 *  We will focus on a cell extending between @f$ x_{\rm cell}-L_x @f$ and @f$ x_{\rm cell} @f$
 *  along the @f$ x @f$-direction, between @f$ y_{\rm cell} - L_y@f$ and @f$ y_{\rm cell} @f$
 *  along the @f$ y @f$-direction, and between @f$ z_{\rm cell} - L_z @f$ and @f$ z_{\rm cell} @f$
 *  along the @f$ z @f$-direction. In the regular grid considered, if @f$(0,0,0)@f$ is taken as 
 *  one of the grid points @f$ x_{\rm cell} @f$, @f$ y_{\rm cell} @f$, and @f$ z_{\rm cell} @f$
 *  will be integer multiples of @f$ L_x @f$, @f$ L_y @f$ and @f$ L_z @f$, respectively.
 * 
 *  In the following we will assume that 
 *  @f[ R_j \textless \min{ \left(L_x, L_y, L_z \right) } @f] 
 *  
 *  The simplest option for relating the positions and masses of the node to 
 *  @f$ \rho_{\rm cell} @f$ is to envision each node @f$ j @f$ as a cube containing @f$ n_j @f$
 *  Kuhn segments, of edge length @f$ R_j @f$, centered at @f$ \mathbf{r}_j @f$.
 *  Node @f$ j @f$ will contribute to a cell if its cube (cube @f$ j @f$) overlaps with the cell. 
 *  Note that, for this to happen, it is not necessary that the nodal position of the center, 
 *  @f$ \mathbf{r}_j @f$, lie in the cell. 
 *  The mass (number of Kuhn segments) contributed by the node to the cell is:
 *  @f[ n_{j,{\rm cell}} = n_j \frac{V_{{\rm cube}\; j\cap{\rm cell}}}{V_{{\rm cube}\;j}} @f]
 *  with @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ being the volume of the intersection of cube 
 *  @f$ j @f$, associated with node @f$ j @f$, and the considered cell, while 
 *  @f$ V_{{\rm cube}\;j} = R_j^3 @f$ is the volume of cube @f$ j @f$.
 *
 * Under the condition @f$ R_j \textless \min{ \left(L_x, L_y, L_z \right) } @f$, 
 *  @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ is obtainable as:
 * \f{eqnarray*}{
 *  V_{{\rm cube}\; j\cap{\rm cell}} & = & 
 *    \max{\left\{\left[ \min{\left(x_j + \frac{R_j}{2}, x_{\rm cell} \right)} 
 *                 -\max{\left(x_j - \frac{R_j}{2}, x_{\rm cell} -L_x \right)}\right] 
 *         , 0 \right\}} \\
 * & \times & 
 *    \max{\left\{\left[ \min{\left(y_j + \frac{R_j}{2}, y_{\rm cell} \right)} 
 *                  -\max{\left(y_j - \frac{R_j}{2}, y_{\rm cell} -L_y \right)}\right] 
 *         , 0 \right\}} \\
 * & \times & 
 *    \max{\left\{\left[ \min{\left(z_j + \frac{R_j}{2}, z_{\rm cell} \right)} 
 *                  -\max{\left(z_j - \frac{R_j}{2}, z_{\rm cell} -L_z \right)}\right] 
 *         , 0 \right\}} 
 * \f}
 * 
 *  As defined by the above equation, @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ is a 
 *  linear function of the node coordinates. Clearly, if cube @f$ j @f$ lies entirely 
 *  within the cell, @f$ V_{{\rm cube}\; j\cap{\rm cell}} = V_{{\rm cube}\:j} @f$ and, 
 *  consequently, @f$ n_{j, {\rm cell}} = n_j @f$. If however, the borders of cube 
 *  @f$ j @f$ intersect the borders of the considered cell, then node @f$ j @f$ will 
 *  contribute a mass @f$ n_{j,{\rm cell}} \textless n_{j} @f$ to the cell. The total mass 
 *  contributed by bead @f$ j @f$ to all cells in which it participates will always 
 *  be @f$n_{j}@f$.
 *
 *  The density @f$ \rho_{\rm cell} @f$ in the considered cell is estimated as: 
 *  @f[ \rho_{\rm cell} = \frac{1}{V_{\rm cell}^{\rm acc}} \sum_j n_{j,{\rm cell}} @f]
 *  Clearly, only nodal points @f$ j @f$ whose cubes have a nonzero overlap with the 
 *  considered cell will contribute to the above summation. The positions vectors 
 *  @f$ \mathbf{r}_j @f$ of these beads will necessarily lie within the considered cell
 *  or its immediate neighbors. 
 * 
 *  The precise conditions for cube @f$ j @f$ to have common points with the considered 
 *  cell are:
 *  @f[ x_{\rm cell} - L_x \textless x_j + \frac{R_j}{2} \textless x_{\rm cell} + R_j @f]
 *  @f[ y_{\rm cell} - L_y \textless y_j + \frac{R_j}{2} \textless y_{\rm cell} + R_j @f]
 *  @f[ z_{\rm cell} - L_z \textless z_j + \frac{R_j}{2} \textless z_{\rm cell} + R_j @f]
 * 
 *  According to the above approach, the force on node @f$ j @f$ due to nonbonded 
 *  interactions is:
 *  @f[ \mathbf{F}_j = - \nabla_{\mathbf{r}_j} A_{\rm nb} 
 *  = - \sum_{\substack{{\rm cells\; having\; common} \\ {\rm points\: with \: cube\;} j}}
 *  V_{\rm cell}^{\rm acc} \left. \frac{d f}{d \rho}\right|_{\rho = \rho{\rm cell}} 
 *  \nabla_{\mathbf{r}_j} \rho_{\rm cell}
 *  @f]
 *    
 *  @section cg_dynamics Coarse-Grained Dynamics
 *  Dynamics at this level of representation will be tracked as two types of processes occuring in 
 *  paralell: (a) Brownian motion of the beads (including crosslinks and end-points) in continuous
 *  three-dimensional space; hops of the slip-springs between adjacent segments along a chain, 
 *  destruction and creation of slip-springs. Both types of processes are governed by the free 
 *  energy of the system.
 * 
 *  @subsection bd_dynamis Brownian Dynamics
 *  We adopt a Brownian Dynamics approach, in which the time evolution of the configuration 
 *  @f$ \left \{ R_j \right\} @f$ is governed by \cite MolPhys_45_637 :
 *  @f[ R_j(t_n + \Delta t) = R_j(t_n) + \frac{\Delta t}{k_{\rm B} T} \: D_{\rm t} \: F_j(t_n) 
 *     + \frac{1}{2} \frac{\Delta t^2}{k_{\rm B} T} \: D_{\rm t} \: \dot{F}_j(t_n)
 *     + X_j^{\rm t}(\Delta t) @f]
 *  where @f$ F_j @f$ is the total force acting on the @f$j@f$-th degree of freedom, 
 *  @f$ \Delta t @f$ the integration timestep and @f$ D_{\rm t} @f$ the translational diffusion 
 *  coefficient of the beads. @f$ X_j^{\rm t} @f$ represents the stochastic 
 *  displacement due to the influence of the coarse-grained microscopic degrees of freedom and 
 *  satisfies the fluctuation-dissipation relation. Random displacements @f$ X_j^{\rm t} @f$ are 
 *  sampled from a Gaussian distribution with zero mean and width: 
 *  @f[ \left\langle (X_j^{\rm t})^2 \right\rangle = 2 D_{\rm t} \Delta t =  
 *     2 \frac{k_{\rm B}T}{m_{\rm bead}\gamma} \Delta t @f]
 *  with @f$ \gamma @f$ being the friction coefficient.
 *  
 *  For large values of @f$\gamma \Delta t @f$ in the diffusive regime, the friction is 
 *  so strong that the velocities relax within @f$\Delta t@f$. For polyisoprene
 *  @f$\zeta_0(500\:{\rm K}) = m_{\rm monomer} \gamma = 2.63 \times 10^{-12} \;{\rm kg/s}@f$, 
 *  @f$\gamma(500\:{\rm K}) = 2.325 \times 10^{13}\; {\rm s}^{-1}@f$ and thus 
 *  @f$\gamma \Delta t=23.25@f$ for @f$ \Delta t = 10^{-12} \;{\rm s}@f$.
 *  
 *  @subsection rotBD Rotational Brownian Dynamics
 *  As originally formulated, the Brownian Dynamics algorithm of Ermak and McGammon 
 *  \cite JChemPhys_69_1352 deals only with the translational motion of the particles. In reality,
 *  however, dispersed particles or molecules also execute rotational Brownian motion arising from 
 *  the fluctuating torque exerted on them by the surrounding solvent molecules. 
 * 
 *  In complete analogy to the Brownian Dynamics equation of motion written for the translational 
 *  degrees of freedom, a similar time evolution equation can be written for the rotational degrees
 *  of freedom  @cite JChemSocFaradayTrans2_81_591 . Thus, angular displacements @f$ \phi_j @f$ are
 *  given by:
 *  @f[ \phi_j\left(t_n+\Delta t\right) = \phi \left(t_n \right) + \frac{\Delta t}{k_{\rm B} T}
 *     \: D_{\rm r} \: T_j \left(t_n \right) + \frac{1}{2} \frac{\Delta t^2}{k_{\rm B} T} \: 
 *     D_{\rm r} \: \dot{T}_j + X_{j}^{\rm r} \left(t_n \right) @f]
 *  where now @f$ D_{\rm r} @f$ stands for the rotational diffusion coefficient of the particles, 
 *  measured in @f$ {\rm rad}^2 \: {\rm s}^{-1} @f$, and @f$ T_{\rm j} @f$ is the sum of external
 *  and interparticle torques acted in direction @f$ j @f$. The timestep @f$ \Delta t @f$ should 
 *  be sufficiently large for angular velocities correlations to vanish out during it.
 * 
 *  @f$ X_j^{\rm r} @f$ represents the stochastic rotation due to the influence of the 
 *  coarse-grained microscopic degrees of freedom and satisfies the fluctuation-dissipation 
 *  relation. Random rotations @f$ X_j^{\rm r} @f$ are 
 *  sampled from a Gaussian distribution with zero mean and width: 
 *  @f[ \left\langle (X_j^{\rm r})^2 \right\rangle = 2 D_{\rm r} \Delta t @f] 
 *  The rotational diffusion coefficient can be related to the rotational frictional coefficient,
 *  @f$ f_{\rm r} @f$, by the Einstein - Smoluchowski equation:
 *  @f[ D_{\rm r} = \frac{k_{\rm B}T}{f_{\rm r}} @f]
 *  The rotational frictional drag coefficient for a sphere of radius @f$ \alpha @f$ is:
 *  @f[ f_{\rm r,sphere} = 8 \pi \eta \alpha^3 @f]
 *  with @f$ \eta @f$ being the dynamic (or shear) viscosity of the medium. 
 *  
 *  @subsection kMC_hopping Kinetic Monte Carlo Hopping
 *  In order to develop a formalism of elementary events of slip-spring hopping, creation or 
 *  destruction, we need expressions for the rate of slippage along the chain backbone. In 
 *  order to extract the diffusivity of the slip-springs, we will proceed along the lines of 
 *  Terzis and Theodorou's work. @cite Macromolecules_35_508 We describe self-diffusion along 
 *  the chain contour with the Rouse model. The Rouse model addresses the dynamics of polymers 
 *  in unentangled melts. A polymer chain is represented by a set of beads connected by 
 *  harmonic springs. The dynamics, as in our simulations, are modeled as a Brownian motion of 
 *  these tethered beads, the environment of a chain being represented as a continuum 
 *  (viscous medium), ignoring all excluded volume and hydrodynamic interactions.
 * 
 *  In this model the self-diffusion of the center of the mass of the polymer is related to 
 *  the friction coefficient, @f$\zeta@f$ on a bead by: 
 *  @f[ D_{\rm Rouse} = \frac{k_{\rm B}T}{N \zeta} @f]
 *  with @f$N@f$ being the number of beads per chain. In the picture we invoke in our network 
 *  model, the center of mass diffusivity along the contour is related to the rate of 
 *  slip-spring jumps across beads (by distance
 *  @f$\left(n_{\rm Kuhns/bead}b^2 \right)^{1/2}@f$ in each direction by (see below for the 
 *  definition of @f$\nu_{\rm diff}@f$)
 *  @f[D_{\rm Rouse} = k_{\rm diff} \frac{n_{\rm Kuhns/bead} b^2}{N} 
 *  = \nu_{\rm diff} \frac{n_{\rm Kuhns/bead} b^2}{N}
 *  \exp{\left(-\frac{A_0}{k_{\rm B}T}\right)} @f]
 *
 *  Hence, one must have: 
 *  @f[\nu_{\rm diff} = \frac{k_{\rm B}T}{n_{\rm Kuhns/bead}b^2 \zeta} 
 *  \exp{\left(-\frac{A_0}{k_{\rm B}T}\right)} @f]
 *  where @f$ A_0 @f$ is a free energy per slip-spring in the equilibrium melt, which 
 *  establishes a baseline for measuring free energies.
 *
 *  An individual jump of the one end of a slip-spring along the chain backbone takes place 
 *  with rate: @f[k_{\rm hopping} = \nu_0 \exp{\left(-\frac{A_{1\to 2}^{\dagger} 
 *  - A_{\rm a_0 - b_0}}{k_{\rm B} T} \right)} 
 *  = \nu_{\rm diff} \exp{\left(-\frac{A_{\rm a_0 - b_0}}{k_{\rm B} T}\right)} @f]
 *  conforming to a transition state theory (TST) picture of the slippage along the backbone 
 *  as an infrequent event, which involves a transition from state "1" to state "2" over a 
 *  free energy barrier. In the final expression, the rate of hopping, 
 *  @f$k_{\rm hopping}@f$, depends directly on the energy of the initial state of the 
 *  slip-spring, @f$A_{\rm a_0 - b_0} @f$, while the dependence on the height of the free 
 *  energy barrier (i.e. @f$A_{1 \to 2}^{\dagger} @f$) has been absorbed into the 
 *  pre-exponential factor, @f$\nu_{\rm diff}@f$.
 * 
 *  @image html hopping.png "Slip-spring hopping schematic"
 *  @image latex hopping.png "Slip-spring hopping schematic" width=15cm
 *
 *  The destruction of a slip-spring can be envisioned as an infrequent event during 
 *  which a slip-spring is pushed to the end of the chain (by traveling a distance 
 *  @f$ \left(n_{\rm Kuhns/bead} b^2 \right)^{1/2}@f$,
 *  corresponding to the last bead of the chain) with rate @f$\nu_{\rm diff}@f$ and 
 *  there faces a transition which can take place with overall rate: 
 *  @f[ k_{\rm destruction} = \nu_{\rm diff} 
 *    \exp{\left(-\frac{A_{\rm a_0-b_0}}{k_{\rm B} T} \right)} @f]
 *  where @f$ A_{\rm a_0 - b_0} @f$ is the free energy the slip-spring (harmonic spring) 
 *  contributes to the total free energy of the system. Again, @f$ \nu_{\rm diff} @f$ is
 *  calculated a la Rouse.
 *
 *  At the timestep of the 3D Brownian Dynamics simulation, when hopping kinetic Monte Carlo has to 
 *  take place, every free end of the system can randomly create a new slip-spring with an internal 
 *  bead of a neighboring chain. This may be accomplished by a rate constant @f$k_{\rm creation} @f$.
 *  The rate for the creation of a slip-spring is closely related to the probability of pairing the 
 *  end ``a'' with one of its candidate mates which lie inside a sphere of prescribed radius 
 *  @f$R_{\rm attempt} @f$. The definition of the probability implies that the more crowded chain ends
 *  are the more probable to create a slip-spring. The number of neighbors around a chain end can be
 *  tuned via the radius of the sphere within which the search takes place, @f$ R_{\rm attempt} @f$. 
 *  A good estimate of @f$R_{\rm attempt} @f$ for polyisoprene (either pure or crosslinked) can be 
 *  given by the tube diameter of the polymer. A computational study of the tube diameter of the 
 *  polyisoprene as a function of the molecular weight has been done by the Li et al.
 *  @cite Polymer_52_5867
 *  The rate constant @f$k_{\rm creation}@f$ can be treated as an adjustable parameter of our model, 
 *  which will be used to ensure that the average number of slip-spring present in the system is 
 *  conserved throughout the simulation.
 * 
 *  @image latex creation_destruction.png "Slip-spring destruction and creation" width=15cm
 * 
 * 
 * @section stress Stress Tensor Calculation
 * Given the free energy functional described in the ``Model'' subsection, the stress tensor of the 
 * system can be derived as \cite Macromolecules_31_6310 :
 * @f[ \bm \sigma = \rho \bm F \left(\frac{\partial A}{\partial \bm F} \right)^{\rm T} @f]
 * where @f$ \bm F @f$ denotes the deformation gradient tensor and may be considered as a mapping of
 * an infinitesimal vector @f$ {\rm d} \bm x @f$ of the initial configuration onto the infinitesimal
 * vector @f$ {\rm d}{\bm x} @f$ of the distorted configuration.
 * 
 * @subsection bond_stress Bonded contributions to the stress tensor
 * The contribution of the bonds to the stress tensor can be easily calcualated due to the fact that
 * invokes central forces between the beads. 
 * The stress tensor of the atom @f$ i @f$, @f$\sigma_{i,\alpha\beta} @f$ is given by the 
 * following formula, where @f$ \alpha @f$ and @f$ \beta @f$ take on values @f$ x @f$, 
 * @f$ y @f$, @f$ z @f$ to generate the six components of the symmetric tensor:
 * @f[ \sigma_{i,\alpha\beta} = -\frac{1}{2}\sum_{j=1}^{N_{\rm b}(i)}
 *      \left(r_{i,\alpha} - r_{j,\alpha}\right)^{\rm min.im.} F_{ij,\beta}^{\rm min.im.}  @f]
 *  where @f$N_{\rm b}(i)@f$ stands for the number of bonds atom @f$ i @f$ participates to and 
 * @f$ F_{ij} @f$ is the force between atoms @f$ i@f$ and @f$ j @f$ calculated by the partial 
 *  derivative of the free energy, @f$ \partial A_{\rm b} / \partial \mathbf{R}_{j} @f$.
 * 
 * @subsection nonbonded_stress Polymer Non-bonded interaction contributions to the stress tensor
 * As already pointed out, the non-bonded intearactions are computed by passing an orthogonal grid
 * in the simulation box and therefore the derivatives present in the definition of the stress
 * tensor are written as a sum over the grid cells:
 * @f[ \frac{\partial A_{\rm nb}}{\partial \bm F} = 
 *     \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(1)}} 
 *     \frac{\partial \rho_{\rm cell}^{(1)}}{\partial \mathbf{F}}
 *    + ... +    
 *     \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(N_{\rm cells})}} 
 *     \frac{\partial \rho_{\rm cell}^{(N_{\rm cells})}}{\partial \mathbf{F}}
 * @f]
 * The derivative of the determinant of @f$ \bm F @f$ with respect to the tensor @f$ \bm F @f$ 
 * itself is calculated by the following equation:
 * @f[ \frac{\partial \left(\det{\left( \mathbf{F} \right)} \right)}{\partial \mathbf{F}}
 *     =  \det{\left(\mathbf{F} \right)} \mathbf{F}^{-1} = \det{\left(\mathbf{F} \right)}
 *     \left( \mathbf{F}^{-1} \right)^{\rm T}
 * @f]
 * 
 * Besides, the determinant of the deformation gradient tensor is written as the ratio of volumes 
 * or densities of the distorted and initial configurations:
 * @f[ \det{\left(\mathbf{F} \right)} = \frac{V^\prime}{V} = \frac{\rho}{\rho^\prime} @f]
 * The symbol ``@f$^\prime @f$'' indicates the initial (undistorted) configuration. Thus, we have:
 * @f[ \frac{\partial \rho_{\rm cell}^{(1)}}{\partial \mathbf{F}} = \rho_{\rm cell}^{\prime \:(1)}
 *  = \rho_{\rm cell}^{\prime \: (1)} \det{\left( \mathbf{F} \right)}
 *  \left(\mathbf{F}^{-1} \right)^{\rm T}
 * @f]
 * and the derivative of the free energy with respect to the local density of a cell becomes:
 * @f[ \left( \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(1)}} 
 *      \frac{\partial \rho_{\rm cell}^{(1)}}{\partial \mathbf{F}} \right)^{\rm T}
 *   = \rho_{\rm cell}^{(1)}\frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(1)}}
 *     \mathbf{F}^{-1}
 * @f]
 * 
 * Similar expressions hold for the rest derivatives and thus, the final form of the stress tensor,
 * due to non-bonded interactions is:
 * @f[ \left(\frac{\partial A_{\rm nb}}{\partial \mathbf{F}} \right)^{\rm T} = 
 *    \left( \rho_{\rm cell}^{(1)} \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(1)}} 
 *    +... + \rho_{\rm cell}^{(N_{\rm cell})} 
 *           \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{N_{\rm cells}}} \right)
 *    \mathbf{F}^{-1} 
 * @f]
 * and
 * @f[ \bm{\sigma}_{\rm nb} = 
 *    \left[ \rho_{\rm cell}^{(1)} \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{(1)}} 
 *    +... + \rho_{\rm cell}^{(N_{\rm cell})} 
 *           \frac{\partial A_{\rm nb}}{\partial \rho_{\rm cell}^{N_{\rm cells}}} \right] 
 *    \mathbf{I}_{3\times 3} @f]
 * It is concluded from the last equation that all off-diagonal elements are equal to zero. Thus, 
 * non-bonded interaction do not contribute to the shear elements of the stress tensor of our model.
 * 
 * @subsection pp_stress Polymer-particle and Particle-particle interaction contributuions to the stress tensor
 * Since polymer-particle and particle-particle interactions yield central forces, the stress tensor
 * contribution can readily be written as:
 * @f[ \sigma_{i,{\rm ps},\alpha\beta} = -\frac{1}{2}\sum_{j=1}^{N_{\rm particles}}
 *      \left(r_{i,\alpha} - r_{j,\alpha}\right)^{\rm min.im.} F_{ij,\beta}^{\rm min.im.}  @f]
 *  where @f$N_{\rm particles}@f$ stands for the number of particles present in the system, and 
 *  @f$ F_{ij} @f$ is the force acted between the @f$ i@f$-th bead and the @f$j @f$-th particle.
 *  The same can also apply for the estimation of the contribution to the stress tensor due to 
 *  particle-particle interactions.
 * 
 * @subsection rheol Estimation of Rheological Properties
 * The main reason for computing the stress tensor is the estimation of rheological properties such
 * as zero-shear viscosity, complex modulus @f$ G^*\left(\omega \right) @f$, storage modulus 
 * @f$ G^\prime \left(\omega\right) @f$ and loss modulus 
 * @f$ G^{\prime \prime} \left(\omega\right) @f$. 
 * Zero-shear viscosity can be calculated as the limit of the off-diagonal stress components 
 * autocorrelation function:
 * @f[ \eta_0 = \frac{V}{k_{\rm B}T} \mathlarger \int_{0}^{+\infty} 
 * \left\langle \sigma_{\alpha \beta}\left(t\right) \sigma_{\alpha \beta}\left(0 \right)
 * \right \rangle {\rm d}t 
 * @f]
 * where @f$ \alpha@f$, @f$ \beta @f$ should be two orthogonal axes.
 * The complex modulus is given as: 
 * @f[ G^* \left(\omega\right) = G^\prime \left(\omega\right) + 
 *    i G^{\prime \prime}\left(\omega\right) = i \omega \frac{V}{k_{\rm B}T}
 *   \mathlarger \int_{0}^{+\infty} {\rm e}^{-i \omega t} 
 * \left\langle \sigma_{\alpha \beta}\left(t\right) \sigma_{\alpha \beta}\left(0 \right)
 * \right \rangle {\rm d}t @f]
 * 
 * @section bench Benchmark Simulations
 * To characterize the equilibrium dynamical behavior, we compute the beads' mean-squared 
 * displacement, @f$ g_1(t) @f$, defined by: 
 * @f[ g_1(t) = \left \langle \left[\mathbf{r}(t) \mathbf{r}(0) \right] \right \rangle @f]
 * where the brackets indicate an average over all beads in the simulation box. We also compute the
 * mean-squared center-of-mass displacements 
 * @f[ g_3(t) = \left \langle \left[ \mathbf{r}_{\rm CM}(t) - \mathbf{r}_{\rm CM}(0) 
 *     \right]^2 \right \rangle @f]
 *
 * We begin our discussion by examining the behavior of a simple melt. To a first approximation, the
 * dynamics of short chains in a melt can be described by the Rouse model @cite JChemPhys_21_1272 . 
 * We have performed simulations for a polymerization index of @f$ N = 80 @f$ beads per chain, each
 * bead representing 10 Kuhn segments of PI (or equivalently 21 PI monomers) having a mass of 
 * @f$m_{\rm bead} = 1460 \:{\rm g/mol} @f$. The simulation box size was fixed at a volume of 
 * @f$ \left(R_{\rm e,0} \right)^3 @f$ with @f$R_{\rm e,0} @f$ being the square root of the
 * unperturbed mean-squared end-to-end distance of the chains and the number of chains present in 
 * the simulation box was set to @f$ n = 27 @f$ in order to match PI's density. The bonded 
 * interactions were parametrized based on the Kuhn length of the PI, @f$ b = 0.9 \;{\rm nm} @f$, 
 * while the non-bonded interactions were parametrized by using the Sanchez - Lacombe equation of \
 * state with PI parameters @cite JSupercritFluids_35_175 .
 *  
 * @subsection rouse Rouse Dynamics
 * The next figure shows results for the time dependence of the mean-squared displacements, 
 * @f$ g_1(t) @f$ and @f$ g_3(t) @f$ for an unentangled melt. 
 * As can be seen, the mean-squared center-of-mass
 * displacement, @f$ g_3(t) @f$, remains linear all times; this means the intermolecular forces 
 * between polymers are too weak to affect diffusive behavior and play a minor role compared to 
 * the bonded interactions. The beads mean-squared displacements, @f$ g_1(t) @f$, exhibit a 
 * subdiffusive behavior that arises from chains' connectivity, and is characterized by a power law
 * of the form @f$ g_1(t) \sim t^{-1/2} @f$. After an initial relaxation time where a change in 
 * @f$ g_1 (t) @f$ occurs, a regular diffusive regime is entered, where @f$ g_1(t) \sim t @f$. This
 * sequence of scaling trends is predicted by the Rouse theory. The limiting behavior of the 
 * chains' center-of-mass displacement can yield an estimate of the diffusivity of the chains: 
 * @f[ \lim_{t\to \infty}{g_3 (t)} = 6 D t @f] 
 * which has been found in excellent agreement with the predicted diffusivity by the Rouse model:
 * @f[ D_{\rm Rouse} = \frac{k_B T}{N \zeta_{\rm bead}} = 1.56 \times 10^{-12} \;
 *     \frac{{\rm m}^2}{s} @f]
 * where our bead's friction coefficient, @f$ \zeta_{\rm bead} @f$ is 
 * @f$ 5.52 \times 10^{-11} \; {\rm kg/s} @f$. The introduction of nonbonded interactions does not
 * seem to affect the scaling laws of the unentangled melt.
 * 
 * @image latex msd_Rouse.eps "Time evolution of the mean-squared displacement of beads and center of mass of the chains for an unentangled PI melt."
 * 
 * @subsection entangled Entangled Dynamics
 * For highly entangled polymer melts, the tube model offers concrete predictions regarding the 
 * scaling behavior of the mean-squared displacement of beads and of the center-of-mass of the 
 * chains. For a very short time the segment does not feel the constraints of the network formed by 
 * the neighboring chains, so that @f$ g_1(t) @f$ is the same as that calcualted for the Rouse model 
 * in free space. Hence, @f$ g_1(t) @f$ can be approximated as:
 * @f[ g_1(t) = \left(\frac{k_{\rm B}T \left(n_{\rm Kuhns/bead} b^2  \right) }
 *     {\zeta_{\rm bead}} \right)^{\frac{1}{2}} \;t^{\frac{1}{2}} @f]
 * This formula is correct when the average displacement is much less than the tube diameter. Let
 * @f$ \tau_{\rm e} @f$ be the time at which the segmental displacement becomes comparable to the 
 * tube diameter @f$ \alpha_{\rm pp} @f$ :
 * @f[ \tau_{\rm e} \simeq \frac{\alpha_{\rm pp}^4 \zeta_{\rm bead}}
 *     {k_{\rm B}T \left(n_{\rm Kuhns/bead} b^2  \right)} @f]
 * The time @f$ \tau_{\rm e} @f$ denotes the onset of the effect of the tube constraints: for 
 * @f$ t \textless \tau_{\rm e} @f$, the chain behaves as a Rouse chain in free space, while for
 * @f$ t \textgreater \tau_{\rm e} @f$ the chain feels the constraints imposed by the tube. 
 * 
 * For @f$ t \textgreater \tau_{\rm e} @f$ the motion of the Rouse segment perpendicular to the 
 * primitive path is restricted, but the motion along the primitive path is free. The mean-squared 
 * displacement along the tube can be approximated as: 
 * @f[ g_1(t) =  
 *    \begin{cases}
 *    \left( \frac{\alpha_{\rm pp}^4  k_{\rm B}T \left(n_{\rm Kuhns/bead} b^2  \right)}
 *    {\zeta_{\rm bead}}\right)^{\frac{1}{4}}  \; t^{\frac{1}{4}}
 *    & , \tau_{\rm e} \lesssim t \lesssim \tau_{\rm R} \\
 *    \left(\frac{ \alpha_{\rm pp} ^2 k_{\rm B}T}{N \zeta_{\rm bead}}  \right)^{\frac{1}{2}}
 *     \; t^{\frac{1}{2}}
 *    &, \tau_{\rm R} \lesssim t \lesssim \tau_{\rm d} 
 *   \end{cases}
 * @f]
 * where the characteristic times are: 
 * @f[ \tau_{\rm R} = \frac{\zeta_{\rm bead} N^2 \left(n_{\rm Kuhns/bead}b^2\right)} 
 *    {3 \pi^2 k_{\rm B}T}
 * @f]
 * and: 
 * @f[ \tau_{\rm d} = \frac{\zeta_{\rm bead} N^3 \left(n_{\rm Kuhns/bead}b^2\right)^2 }
 *     {\pi^2 k_{\rm B}T \alpha_{\rm pp}^2} 
 * @f]
 * 
 * For @f$ t \textgreater \tau_{\rm d} @f$, the dynamics is governed by the reptation process. 
 * The mean-squared displacement of a bead is approximated by:
 * @f[ g_1(t) \simeq \frac{k_{\rm B}T\alpha_{\rm pp}^2} {N^2 \zeta_{\rm bead} 
 *   \left(n_{\rm Kuhns/bead} b^2 \right)} \: t 
 * @f]
 *  
 * In order to reproduce the entangled dynamics of the PI melts we start from the unentangled 
 * melt of the previous section and we gradually introduce slip-springs between the chains by 
 * following the kinetic Monte Carlo algorithm described. We can tune the rate of slip-spring 
 * creation process, @f$ k_{\rm creation} @f$ in a way that the system stabilizes at the experimental
 * count of entanglements, as expected by the average entanglement molecular weight. 
 * This procedure is depicted in the following figure, where our simulation
 * starts with a completely unentangled melt until it reaches a fully entangled system, 
 * representative of the experimental PI.
 * 
 * @image latex sl_count.eps "Number of slip-springs present in the system as a function of time."
 * 
 * Finally, we study the mean-squared displacements @f$ g_1(t) @f$ and @f$ g_3(t) @f$ as a function 
 * of time, by using the simulation trajectory after the point where the number of slip-sprins
 * has stabilized at its mean value. It can be seen in the following figure that at short
 * times the bead's mean-squared displacement, @f$ g_1(t) @f$, shows a scaling regime with a power
 * law @f$ t^{1/2} @f$; at intermediate times, a regime with a power law @f$ t^{1/4} @f$ appears and
 * eventually we observe a crossover to regular diffusion at long times. Before the diffusive 
 * behavior appears, tube model predicts a crossover to @f$ t^{1/2} @f$ which is completely 
 * absent from our model. Together with the simulation results, the experimental estimations of the
 * the characteristic times of PI are presented in the figure. Our model seems capable of boldly
 * reproducing the dynamics of entangled melts. However, careful parametrization is needed before
 * rheological predictions can be extracted.
 * The mean-squared displacement of the chains center-of-mass, @f$ g_3(t) @f$, also exhibits 
 * subdiffusive behavior 
 * at intermediate times, with a scaling behavior of @f$ t^{1/2} @f$, as predicted by the tube 
 * model; at long times regular diffusion is achieved. 
 * 
 * @image latex msd_gaussian_entangled.eps "Time evolution of the mean-squared displacement of beads and center of mass of the chains for an entangled PI melt."
 * 
 */
