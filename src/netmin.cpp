/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @version 1.1 (November 28, 2013)
 *  @brief The routines of the "NetwMin" application class are defined here.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <string>
#include "netmin.h"

using namespace std;

namespace NetworkNS{
   
   /// @param[in] filename The name of the data file to be read in order to initialize the 
   ///                     simulation.
   NetwMin::NetwMin(std::string filename){
      
      domain = new Domain(filename);
      cout << "#: Domain characteristics have been read.\n";
      network = new Network(this,filename);
      cout << "#: Network has been read.\n";
      my_rnd_gen = new RanMars(11111);
      
      cout << "#: Particles have been read.\n";
      
      my_traj_file = new dump("dump_b3D.lammpstrj");
      cout << "#: 'dump_b3D.lammpstrj' will be used as trajectory file.\n";
   
      grid = new Grid(domain->XBoxLen, domain->YBoxLen, domain->ZBoxLen, 5, 5, 5);
      write_network_to_lammps_data_file();
      
      return;
   }
   

   void NetwMin::write_network_to_lammps_data_file() {

      double wrapped_coords[3];
      int  pcoeffs[3];
      
      pcoeffs[0] = 0;
      pcoeffs[1] = 0;
      pcoeffs[2] = 0;
      
      FILE *lammps_file = fopen("restart.data", "wt");


      /* Write header to LAMMPS data file: */
      
      fprintf(lammps_file, "LAMMPS data file of crosslinked coarse-grained network.\n\n");
      fprintf(lammps_file, "%ld atoms\n", network->nodes.size());
      fprintf(lammps_file, "%ld bonds\n", network->strands.size());
      fprintf(lammps_file, "\n");
               
      fprintf(lammps_file, "%lu atom types\n", network->node_types.size());
      
      if (network->pslip_springs.size() > 0) 
         fprintf(lammps_file, "2 bond types\n");
      else
         fprintf(lammps_file, "1 bond types\n");
      fprintf(lammps_file, "\n");

      fprintf(lammps_file, "%lf %lf xlo xhi\n", 
              domain->BoxLow[0], domain->BoxHigh[0]);
      fprintf(lammps_file, "%lf %lf ylo yhi\n", 
              domain->BoxLow[1], domain->BoxHigh[1]);
      fprintf(lammps_file, "%lf %lf zlo zhi\n", 
              domain->BoxLow[2], domain->BoxHigh[2]);
      
      
      fprintf(lammps_file, "\nMasses\n\n");
      for (unsigned int i = 0; i < network->node_types.size(); i++)
         fprintf(lammps_file, "%d \t %lf \t %-1.4e \t # n_Kuhns/bead - g/mol\n", 
                 i+1, network->node_types[i].n_mass, network->node_types[i].mass);
      
      
      fprintf(lammps_file, "\nBond Coeffs\n\n");
      for (unsigned int i = 0; i < network->bond_types.size(); i++)
         fprintf(lammps_file, 
                 "%d \t %-1.4e \t %-1.4e \t %-1.4e \t # 3/2*k_B (kJ/mol/K) - <R_e^2> (A^2) - b (A)\n",
                 i+1, network->bond_types[i].spring_coeff, 
                 network->bond_types[i].sq_ete, network->bond_types[i].kuhnl);
      
      
      fprintf(lammps_file, "\nPair Coeffs\n\n");
      for (unsigned int i = 0; i < network->node_types.size(); i++)
         fprintf(lammps_file, "%d \t %lf \t # r_node = b*n_j^(1/2)\n",i+1,
                 network->node_types[i].r_node);

      
      fprintf(lammps_file, "\nAtoms\n\n");
      int inode = 0;
      for (std::list<tNode> ::iterator it = network->nodes.begin(); 
           it != network->nodes.end(); ++it) {
         inode++;
         (*it).Id = inode;
         
         wrapped_coords[0] = (*it).Pos[0];
         wrapped_coords[1] = (*it).Pos[1];
         wrapped_coords[2] = (*it).Pos[2];
         //domain->put_in_primary_box(wrapped_coords, pcoeffs);
         
         fprintf(lammps_file, "%d  %d  %d  0.000  %.10f  %.10f  %.10f  %d  %d  %d\n",
                 (*it).Id, (*it).OrChains[0] ,(*it).Type, 
                 wrapped_coords[0], wrapped_coords[1], wrapped_coords[2],
                 pcoeffs[0], pcoeffs[1], pcoeffs[2]);
      }

      fprintf(lammps_file, "\nBonds\n\n");
      
      int istrand = 0;
      for (std::list<tStrand> ::iterator it = network->strands.begin(); 
           it != network->strands.end(); ++it) { 
         istrand++;
         
         if ((*it).slip_spring)
            fprintf(lammps_file, "%d  2  %d  %d # %d\n",
                 istrand, (*it).pEnds[0]->Id, (*it).pEnds[1]->Id, (*it).OrChain);
         else
            fprintf(lammps_file, "%d  1  %d  %d # %d\n",
                 istrand, (*it).pEnds[0]->Id, (*it).pEnds[1]->Id, (*it).OrChain);
      }
      fclose(lammps_file);

      return;
   }

}