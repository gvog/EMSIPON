/**
 * @file
 * @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 * @brief  The class of the network application itself.
 * @version 4.0 (November 28, 2013)
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef NETWMIN_H
#define	NETWMIN_H

#include <vector>

#include "domain.h"
#include "dump.h"
#include "grid.h"
#include "hopping.h"
#include "rng.h"
#include "network.h"


namespace NetworkNS {

   
   /** @class NetwMin
    *  @brief The class of the host application itself. It contains pointers to all constituents, e.g. 
    *         the simulation domain, random number generator, etc. 
    */
   
   class NetwMin {
   public:
      class Network *network;    ///< A pointer to a network class hosting the connectivity of 
                                 ///< the system.
      
      class Domain *domain;      ///< A pointer to a domain class hosting the dimensions of the 
                                 ///< simulation domain.
      
      class Grid *grid;          ///< A pointer to a grid class hosting the information about the 
                                 ///< nonbonded free energy estimation grid.
      
      class RanMars *my_rnd_gen; ///< A pointer to the class of the random number generator.
      
      class dump *my_traj_file;  ///< A pointer to the class which takes care of writing the 
                                 ///< trajectory file of the simulation.
      
      NetwMin(std::string);      ///< The constructor of the application. The only input is the 
                                 ///< name of the data file to read the initial configuration from.
      
      void write_network_to_lammps_data_file(); ///< A function for writing the current 
                                                ///< configuration of the system in LAMMPS-like 
                                                ///< format. It can be used for restarting a 
                                                ///< simulation.
      
   };
}

#endif	/* NETWMIN_H */

