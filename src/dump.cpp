/** 
 * @file
 * @author  Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 * @version 1.0 (created on October 28, 2013)
 * @brief   This file contains all routines necessary to write a trajectory file of the simulation.
 *          It uses a pretty standard LAMMPS trajectory format.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include "dump.h"
#include "network.h"
#include "netmin.h"
#include "b3D_integrator.h"

using namespace std;

namespace NetworkNS {

   dump::dump(std::string filename_to_open) {
      
      my_file.open(filename_to_open.c_str(), ofstream::out);
      if (not my_file.good())
         cout << "Specified dump file does not exist!" << endl;

   }

   /** @param[in] netw_app A pointer to a network application, which contains all relevant 
    *                      information concerning
    *                      the topology of the system under investigation. 
    *  @param[in] b3D      A pointer to a Brownian Dynamics integrator containing the current 
    *                      positions of the bead in the course of the simulation.
    *  @param[in] timestep The current timestep of the simulation. 
    */
   void dump::add_snapshot_to_dump(const class NetwMin *netw_app, const class cb3D_integrator *b3D, 
                                   unsigned int timestep){
      
      my_file << "ITEM: TIMESTEP\n";
      my_file << timestep << "\n";
      my_file << "ITEM: NUMBER OF ATOMS\n";
      my_file << netw_app->network->nodes.size() << "\n";
      my_file << "ITEM: BOX BOUNDS pp pp pp\n";
      my_file << netw_app->domain->BoxLow[0] << " " << netw_app->domain->BoxHigh[0] << "\n";
      my_file << netw_app->domain->BoxLow[1] << " " << netw_app->domain->BoxHigh[1] << "\n";
      my_file << netw_app->domain->BoxLow[2] << " " << netw_app->domain->BoxHigh[2] << "\n";
      my_file << "ITEM: ATOMS id x y z ix iy iz Vor VorNeigh xx yy zz xy xz yz\n";
      
      // TODO: Here we can ask for a Voronoi tessellation of the simulation box in order to 
      //       calculate atomic volumes.
      
      my_file.precision(6);
      for (std::list <tNode>::iterator it = netw_app->network->nodes.begin();
           it != netw_app->network->nodes.end(); ++it){
         
         unsigned int ibead = (*it).Id - 1;
         
         my_file << (*it).Id << " " 
                 << b3D->bd_x[3*ibead + 0] << " "
                 << b3D->bd_x[3*ibead + 1] << " "
                 << b3D->bd_x[3*ibead + 2] << " "
                 << " 0 0 0 1.0 0 " 
                 << b3D->bd_stress[ibead][0] << " "
                 << b3D->bd_stress[ibead][1] << " "
                 << b3D->bd_stress[ibead][2] << " "
                 << b3D->bd_stress[ibead][3] << " "
                 << b3D->bd_stress[ibead][4] << " "
                 << b3D->bd_stress[ibead][5] << endl;
      }
   }
   
   dump::dump(const dump& orig) {
      
   }

   dump::~dump() {
      my_file.close();
   }

}