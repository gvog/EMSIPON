/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @brief  The header file of the slip-spring hopping kinetic Monte Carlo scheme.
 *  @version 2.0 (Created on May 14, 2013)
 */

#ifndef HOPPING_H
#define	HOPPING_H

#include <fstream>
#include <iostream>
#include <string>

#include "b3D_integrator.h"
#include "net_types.h"
#include "network.h"

using namespace std;
using namespace NetworkNS;

namespace NetworkNS{
   
   /** @class Hopping
    *  @brief The class of the hopping kinetic Monte Carlo scheme. It can be called by a Brownian Dynamics
    *         class, get all the necessary information from it and alter the connectivity of the system,
    *         based on the rates described in hopping.cpp .
    */
   class Hopping {
      
   public:
      Hopping(double);        ///< The constructor of the class, where a file for writing the lifetime of the 
                        ///< slipsprings is created. 
      ~Hopping();       ///< The destructor of the class.
      
      void hopping_step(class NetwMin *netapp, const class cb3D_integrator *b3D, double *pos_array, 
                        double temperature, double elapsed_time); 
      ///< The routine which executes a single kinetic Monte Carlo slipspring hopping step.
   
   private:
      std::ofstream lifetimes_file; ///< A file for writing the lifetimes of the slipsprings to.
      
      std::ofstream events_file; ///< A file for writing the hopping events taking place.
      
      double nu_hopping_times_exp_of_barrier;
      
   };
   
   
}




#endif	/* HOPPING_H */

