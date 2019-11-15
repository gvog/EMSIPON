/** @file
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @brief  C++ source file containing the implementation of bonded interactions.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Grigorios Megariotis and Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <cmath>
#include <iostream>

#include "distributions.h"

using namespace std;

/** @param[in] rij the separation vector between beads i and j, i.e. @f$ \mathbf{R}_{ij} @f$.
 *  @param[in] coeff the strength of the entropic springs 
 *                   @f$ \epsilon_{\rm b} = 3/2 \:k_{\rm B} @f$, measured in kJ/mol/K.
 *  @param[in] sq_ete the equilibrium mean-squared end-to-end length of the strand, 
 *                    i.e. @f$\sigma_{\rm b} =  n_{\rm Kuhns/bead}b^2@f$, measured 
 *                    in @f$ \text{\AA}^2 @f$.
 *  @param[in] temp the temperature of the simulation, @f$ T@f$, in K.
 *  @param[out] gradi the force acted on the bead i, due to its bond with bead j
 *  @param[out] gradj the force acted on the bead j, due to its bond with bead i
 */ 
double f_gaussian(const double *rij, const double & coeff, const double & sq_ete, 
              const double & temp, double *gradi, double *gradj) {

   /** This routine applies a Gaussian free energy potential of the form:
    *  @f[ \mathcal{V}_{\rm b} \left(r^2_{ij} \right) = \epsilon_{\rm b} T 
    *  \frac{r^2_{ij}}{\sigma^2_{\rm b}}
    *  = \frac{3}{2}k_{\rm B}T\frac{\mathbf{R}_{ij} \cdot \mathbf{R}_{ij}}
    *  {n_{\rm Kuhns/bead} b^2} @f]
    * with the parameters @f$ \epsilon_{\rm b} @f$ and @f$ \sigma_{\rm b} @f$ read from the 
    * data file. In our approach @f$ \epsilon_{\rm b} = 0.012471 @f$  kJ/mol/K and 
    * @f$ \sigma_{\rm b} = 810 \; \text{\AA}^2 @f$ for polyisoprene melt.
    */
   
   // Compute the distance between positional vectors ri and rj
   double rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

   // Compute the free energy of the Gaussian spring.
   double p = coeff * temp / sq_ete;

   // Compute the forces as gradients of the Gaussian distribution along ri and rj direction.
   gradi[0] = 2.0 * p * rij[0];
   gradi[1] = 2.0 * p * rij[1];
   gradi[2] = 2.0 * p * rij[2];

   gradj[0] = -gradi[0];
   gradj[1] = -gradi[1];
   gradj[2] = -gradi[2];
      
   return (p*rsq);
}


/** @param[in] rij the separation vector between beads i and j, i.e. @f$ \mathbf{R}_{ij} @f$.
 *  @param[in] coeff the strength of the entropic springs 
 *                   @f$ \epsilon_{\rm b} = 3/2 \:k_{\rm B} @f$, measured in kJ/mol/K.
 *  @param[in] sq_ete the equilibrium mean-squared end-to-end length of the strand, 
 *                    i.e. @f$\sigma_{\rm b} =  n_{\rm Kuhns/bead}b^2@f$, measured in 
 *                    @f$ \text{\AA}^2 @f$.
 *  @param[in] temp the temperature of the simulation, @f$ T@f$, in K.
 */ 
double e_gaussian(const double *rij, const double & coeff, const double & sq_ete, 
                  const double & temp) {
   //compute the distance between positional vectors ri and rj
   double rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

   // Compute only the free energy of the Gaussian spring.
   return (coeff * temp * rsq / sq_ete);
}

