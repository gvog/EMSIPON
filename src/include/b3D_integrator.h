/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @version 1.0 (January 15, 2013)
 *  @brief   Header file accompanying the ``b3D_integrator.cpp'' C++ source file.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis and Grigorios Megariotis. 
 *  All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef _B3D_INTEGRATOR_H
#define _B3D_INTEGRATOR_H

#include <cmath>
#include <cstdio>
#include <list>
#include <stdlib.h>
#include <vector>

#include "hopping.h"
#include "net_types.h"
#include "network.h"
#include "domain.h"
#include "rng.h"

namespace NetworkNS {
   
   /** @class cb3D_integrator
    *  @brief The class of the Brownian Dynamics integrator.
    */
   
   class cb3D_integrator{
   public:
      
      cb3D_integrator(class NetwMin *, double, double); ///< The constructor of the class.
      ~cb3D_integrator();                       ///< The destructor of the class.
      
      void integrate(unsigned int nsteps, double dt, unsigned int nstout);  
      ///< The main function which performs the time integration, following Brownian Dynamics.
  
      void extract_positions(double *x); ///< A function for extracting the positions from the integrator 
      
      void compute_stresses (void);    ///< A function for computing the per-atom stress.
      
      void report(unsigned int, double, double);   ///< The function taking care of reporting the current status of 
                                                   ///< the simulation to both the standard output and the dump file.

      double bonded_force_calculation(bool); ///< The function for the calculation of bonded interactions
      
      double simpler_scheme_non_bonded_force_calculation(void);   
      
      void calculate_pressure(double *);  ///< Sum the per-atom stresses to calculate the pressure of the 
                                          ///< simulation box.
      
      double *bd_gamma;    ///< The friction coefficient, @f$\gamma@f$, of every degree of freedom, 
                           ///< measured in @f$ {\rm s}^{-1} @f$.
      
      double *bd_mass;     ///< The mass, @f$m_i @f$, of every degree of freedom, measured in 
                           ///< @f$ {\rm g}/{\rm mol} @f$.
      
      double *bd_x;        ///< The current position of the beads during the Brownian Dynamics integration, 
                           ///< at the current timestep, i.e. @f$ \mathbf{r}_i \left( t\right) @f$. The 
                           ///< size of the vector is three times the degrees of freedom of the BD 
                           ///< simulation.
      double **bd_stress;  ///< The per-atom stress tensor.                      
  
      unsigned int bd_cur_step;  ///< The current step of the integration.
      
   private:
      
      class NetwMin *cur_bd_net; ///< The network application that created the class.

      unsigned int dofs;      ///< The count of particles simulated using the BD scheme.
      unsigned int dofs_3N;   ///< The number of degrees of freedom for which BD takes place.
      
      
    
      double *bd_f;     ///< The forces due to bonded interactions during the Brownian Dynamics 
                        ///< simulation, at the current timestep, i.e. 
                        ///< @f$ \mathbf{F}_i \left(t \right)@f$.
      
      double *bd_nb_f;  ///< The forces due to nonbonded interactions duringn the BD simulations, at the
                        ///< current timestep. The size of the vector is three times the degrees of 
                        ///< freedom of the BD simulation.
      
      double *bd_x_ps;  ///< The positions of the beads of the BD simulation at the previous timestep.
                        ///< The size of the vector is 
      
      double *bd_f_ps;  ///< The forces acted on the beads durign the previous timestep of the BD 
                        ///< simulation.
      
      class Hopping *my_hopping_scheme;   ///< A pointer to a hopping kinetic Monte Carlo class. 
      
      clock_t tbegin;   ///< The intial timestep of the Brownian Dynamics simulation.
      
      double *xshift;   ///< The vector for keeping the x coordinates of the beads with respect to the nonbonded
                        ///< free energy estimation grid, xshift @f$\in [0, L_{\rm x}) @f$ with @f$ L_{\rm x} @f$ 
                        ///< being the @f$x@f$ edge length of the simulation box.
      
      double *yshift;   ///< The vector for keeping the y coordinates of the beads with respect to the nonbonded
                        ///< free energy estimation grid, yshift @f$\in [0, L_{\rm y}) @f$ with @f$ L_{\rm y} @f$ 
                        ///< being the @f$ y @f$ edge length of the simulation box.
      
      double *zshift;   ///< The vector for keeping the z coordinates of the beads with respect to the nonbonded
                        ///< free energy estimation grid, zshift @f$\in [0, L_{\rm z}) @f$ with @f$ L_{\rm z} @f$ 
                        ///< being the @f$z @f$ edge length of the simulation box.
      
      int *grid_cell;   ///< The global tag of the cell the node belongs to.
       
      double bd_temp;   ///< The temperature of the Brownian Dynamics simulation.
          
      bool gamma_mass_opt;    ///< A boolean variable indicating whether or no the BD integrator run in an optimized
                              ///< way. If true, all beads have the same mass and friction coefficient, so the 
                              ///< integrator does not use arrays in order to speed-up the execution.
     
      double *density_cells;    ///< The local density per cell. 
      
      double *den_dx;   ///< The partial derivative of the local cell density with respect to the @f$x@f$-coordinate
                        ///< of a bead, @f$\nabla_{\mathbf{r}_j} \rho_{\rm cell} @f$.
                           
      double *den_dy;   ///< The partial derivative of the local cell density with respect to the @f$y@f$-coordinate
                        ///< of a bead, @f$\nabla_{\mathbf{r}_j} \rho_{\rm cell} @f$.
      
      double *den_dz;   ///< The partial derivative of the local cell density with respect to the @f$z@f$-coordinate
                        ///< of a bead, @f$\nabla_{\mathbf{r}_j} \rho_{\rm cell} @f$.
      
  
      FILE * p_cell_density;
      
      void from_bd_x_to_polymer_network(void);     // Convert elements of array bd_x to positions of the polymer network
     
      void cell_density_nodal_points(void);        // compute the density in each cell of the orthogonal grid    
      
   };
   
}  

#endif	/* NETWORK_H */




