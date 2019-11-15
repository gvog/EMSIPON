/** 
 *  @file 
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @brief The C++ source file containing all functions relevant to the class Network.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Georgios Vogiatzis. 
 *  All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <vector>
#include <cmath>

#include "constants.h"
#include "netmin.h"
#include "Auxiliary.h"


using namespace std;


namespace NetworkNS {

   Network::~Network() {
      nodes.clear();
      strands.clear();
   }

   
   /** @param[in] netw_min A pointer to the application which initializes the constructor. 
    *                      The application is a class NetwMin. 
    *  @param[in] filename A standard C++ string containing the name of the data file to open in 
    *                      order to read the polymeric network from.
    */
   Network::Network(class NetwMin *netw_min, std::string filename) {

      int nnodes = 0, nstrands = 0, atoms_line_start = 0, bonds_line_start = 0,
              bond_coeffs_line_start = 0, masses_line_start = 0, pair_coeffs_line_start = 0,
              atom_types_line = 0, bond_types_line = 0, ihelp;

      nstrands = 0;
      cout << "#: Initializing network.\n";
      /* Define and open the desired file. */
      ifstream data_file(filename.c_str(), ifstream::in);
      /* Define an array of strings to hold the contents of the file. */
      std::vector <string> lines_of_file;
      std::vector <string> tokens;

      /* Read and search file. */
      int iline = 0;

      while (data_file.good()) {
         /* A temporary string for the current line of the file. */
         std::string current_line;
         getline(data_file, current_line);
         /* Add current line to file's array of lines. */
         lines_of_file.push_back(current_line);

         if (current_line.find("atoms") != string::npos)
            nnodes = atoi(current_line.c_str());
         if (current_line.find("bonds") != string::npos)
            nstrands = atoi(current_line.c_str());
         if (current_line.find("Atoms") != string::npos)
            atoms_line_start = iline;
         if (current_line.find("Bonds") != string::npos)
            bonds_line_start = iline;
         if (current_line.find("Masses") != string::npos)
            masses_line_start = iline;
         if ((current_line.find("Bond") != string::npos)
                 && (current_line.find("Coeffs") != string::npos))
            bond_coeffs_line_start = iline;
         if ((current_line.find("Pair") != string::npos)
                 && (current_line.find("Coeffs") != string::npos))
            pair_coeffs_line_start = iline;
         if ((current_line.find("atom") != string::npos)
                 && (current_line.find("types") != string::npos))
            atom_types_line = iline;
         if ((current_line.find("bond") != string::npos)
                 && (current_line.find("types") != string::npos))
            bond_types_line = iline;

         iline++;
      }
      data_file.close();

      /* Read the masses of the nodes from the corresponding session of the data file. */
      tokens = tokenize(lines_of_file[atom_types_line]);
      ihelp = atoi(tokens[0].c_str());
      cout << "#: --- Number of bead types: " << ihelp << " .\n";
      node_types.resize(ihelp);
      for (int i = 0; i < ihelp; i++) {
         tokens = tokenize(lines_of_file[masses_line_start + i + 2]);
         node_types[i].n_mass = atof(tokens[1].c_str());
         node_types[i].mass = atof(tokens[2].c_str());

         tokens = tokenize(lines_of_file[pair_coeffs_line_start+i+2]);
         node_types[i].r_node = atof(tokens[1].c_str());
      }

      /* Read the bond coefficients from the corresponding section of the data file. */
      tokens = tokenize(lines_of_file[bond_types_line]);
      ihelp = atoi(tokens[0].c_str());
      cout << "#: --- Number of bond types: " << ihelp << " .\n";
      bond_types.resize(ihelp);
      for (int i = 0; i < ihelp; i++){
         tokens = tokenize(lines_of_file[bond_coeffs_line_start + i + 2]);
         bond_types[i].spring_coeff = atof(tokens[1].c_str());
         bond_types[i].sq_ete = atof(tokens[2].c_str());
         bond_types[i].kuhnl = atof(tokens[3].c_str());
      }

      /* Read atom sections which corresponds to nodal points: */
      for (int i = 0; i < nnodes; i++) {
         /* Read and split the line from the input file: */
         tokens = tokenize(lines_of_file[atoms_line_start + 2 + i]);

         
         //int pbcx = atoi(tokens[7].c_str());
         //int pbcy = atoi(tokens[8].c_str());
         //int pbcz = atoi(tokens[9].c_str());
         
         /* Set the attributes of the current node object. */
         tNode current_node;
         current_node.Id = atoi(tokens[0].c_str());
         current_node.Type = atoi(tokens[2].c_str());
         current_node.Pos[0] = atof(tokens[4].c_str()); //+ (double)pbcx*netw_min->domain->XBoxLen;
         current_node.Pos[1] = atof(tokens[5].c_str()); //+ (double)pbcy*netw_min->domain->YBoxLen;
         current_node.Pos[2] = atof(tokens[6].c_str()); //+ (double)pbcz*netw_min->domain->ZBoxLen;
         current_node.node_cell = 0;

         // Read the mass of the nodal point from the "Masses" section of the file.
         current_node.n_mass = node_types[current_node.Type-1].n_mass;
         current_node.mass = node_types[current_node.Type-1].mass;

         /* Pair coeffs come from the data file in Angstrom. */
         current_node.r_node = node_types[current_node.Type-1].r_node;
         current_node.r_star = current_node.r_node;

         nodes.push_back(current_node);
      }
      cout << "#: --- " << nnodes << " nodes have been added to the network structure.\n";


      int nchains = 0;

      for (int i = 0; i < nstrands; i++) {
         /* Read and split the line from the input file: */
         tokens = tokenize(lines_of_file[bonds_line_start + 2 + i]);
         
         /* Create a temporary object: */
         tStrand current_strand;
         current_strand.Id = atoi(tokens[0].c_str());
         current_strand.Type = atoi(tokens[1].c_str());
         current_strand.OrChain = atoi(tokens[5].c_str());
         if (current_strand.OrChain > nchains)
            nchains = current_strand.OrChain;

         if (current_strand.OrChain == 0)
            current_strand.slip_spring = true;
         else
            current_strand.slip_spring = false;

         /* Read the tags of the start and the end of the strand. */
         int snode = atoi(tokens[2].c_str());
         int enode = atoi(tokens[3].c_str());

         if (snode > (nnodes) || enode > (nnodes)){
            cout << ": --- There is a problem at line: " << (bonds_line_start + 2 + i + 1) 
                 << " of the data file.\n" << lines_of_file[bonds_line_start + 2 + i] << "\n";
            exit(0);
         }
         
         /* We have to find the nodal point whose tag is equal to snode: */
         int found = 0;
         for (std::list <tNode>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
            if (found > 1)
               break;
            if ((*it).Id == snode || (*it).Id == enode) {
               found++;
               (*it).OrChains.push_back(current_strand.OrChain);
               current_strand.pEnds.push_back(&(*it));
            }
         }

         /* Spring coeffs comes from the data file in kJ/mol/K/A^2 */
         current_strand.spring_coeff = bond_types[current_strand.Type-1].spring_coeff;
         /* Length of the strand (n*b) comes from the data file in Angstroms. */
         current_strand.sq_end_to_end = bond_types[current_strand.Type-1].sq_ete;
         current_strand.kuhn_length = bond_types[current_strand.Type-1].kuhnl;

         /* Set the creation time of the slip-spring: */
         if (current_strand.slip_spring)
            current_strand.tcreation = 0;
         
         /* and push it back to the list: */
         strands.push_back(current_strand);

         /* 2013/04/09: Also, push back a pointer to the last element of "strands" vector, if 
          *             the inserted strand is a slip spring. */
         if (current_strand.slip_spring)
            pslip_springs.push_back(&(strands.back()));
      }



      for (std::list<tStrand> ::iterator it = strands.begin(); it != strands.end(); ++it) {
         /* 2013/04/09: We exclude slip strings from the local registry of strands connected to 
          *             a specific node. Array "pStrands" contains only chain strands, 
          *             not slip-springs!!. */
         if (!(*it).slip_spring) {
            (*it).pEnds[0]->pStrands.push_back(&(*it));
            (*it).pEnds[1]->pStrands.push_back(&(*it));
         }
      }

      for (std::list<tStrand> ::iterator it = strands.begin(); it != strands.end(); ++it) {
         if ((*it).slip_spring)
            continue;

         double dxr = (*it).pEnds[1]->Pos[0] - (*it).pEnds[0]->Pos[0];
         double dyr = (*it).pEnds[1]->Pos[1] - (*it).pEnds[0]->Pos[1];
         double dzr = (*it).pEnds[1]->Pos[2] - (*it).pEnds[0]->Pos[2];
         netw_min->domain->minimum_image(dxr, dyr, dzr);

         double dist = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
         if (dist - (*it).sq_end_to_end / (*it).kuhn_length > tol)
            cout << "Strand " << (*it).Id << " has length " << (*it).sq_end_to_end 
                  / (*it).kuhn_length << " but real dist " << dist << endl;

      }

      /* Categorize strands into chains: */
      std::vector<std::list<tStrand *> > chains(nchains);
      for (std::list<tStrand> ::iterator it = strands.begin(); it != strands.end(); ++it)
         if (!(*it).slip_spring)
            chains[(*it).OrChain - 1].push_back(&(*it));


      /* Loop over all chains: */
      sorted_chains.resize(nchains);

      for (unsigned int ich = 0; ich < chains.size(); ich++) {
         /* Search to find an end of the chain: */
         std::list <tStrand *>::iterator it, jt;
         jt = find_if(chains[ich].begin(), chains[ich].end(), pred_strand_is_end);
         if (jt == chains[ich].end()) {
            cout << "Chain with no ends was found: " << ich << " .\n";
            cout << "Number of strands: " << chains[ich].size() << " .\n";
            return;
         }

         if ((*jt)->pEnds[0]->Type != 1) {
            tNode *temp;
            temp = (*jt)->pEnds[0];
            (*jt)->pEnds[0] = (*jt)->pEnds[1];
            (*jt)->pEnds[1] = temp;
         }
         sorted_chains[ich].push_back((*jt));
         if (chains[ich].size() > 1)
            chains[ich].erase(jt);
         else
            continue;

         int nelems = chains[ich].size();
         for (int i = 0; i < nelems; i++) {
            /* Set the iterator to the end of the sorted list: */
            it = sorted_chains[ich].end();
            --it;

            /* Find the next element: */
            for (jt = chains[ich].begin(); jt != chains[ich].end(); ++jt) {
               /* check whether is connected to the preceeding strand: */
               if ((*jt)->pEnds[0] == (*it)->pEnds[1]) {
                  /* We can append the found strand as is: */
                  sorted_chains[ich].push_back((*jt));
                  /* an delete it from the previous list */
                  jt = chains[ich].erase(jt);
                  break;
               }

               if ((*jt)->pEnds[1] == (*it)->pEnds[1]) {
                  /* We have to change the order of the pointers. */
                  tNode *temp;
                  temp = (*jt)->pEnds[0];
                  (*jt)->pEnds[0] = (*jt)->pEnds[1];
                  (*jt)->pEnds[1] = temp;

                  /* And add it to the list: */
                  sorted_chains[ich].push_back((*jt));
                  jt = chains[ich].erase(jt);
                  break;
               }
            }
         }
      }


      return;
   }


   /// @param[in] current A pointer to a type tStrand. The strand to check whether is located at a 
   ///                    chain end or no.
   bool pred_strand_is_end(tStrand * current) {
      bool outcome;
      if (current->pEnds[0]->Type == 1 || current->pEnds[1]->Type == 1)
         outcome = true;
      else
         outcome = false;

      return outcome;
   }

}
