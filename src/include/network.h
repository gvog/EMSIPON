/**
 * @file
 * @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 * @version 2.1 (May 12, 2012)
 * @brief The header file of the Network class itself, which describes the topology of the system under consideration.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef _NETWORK_H
#define _NETWORK_H

#include "net_types.h"

namespace NetworkNS {

   /** @class Network 
    *  @brief The class which stores all information concerning the polymeric network.
    */
   class Network {
   public:
      
      Network(class NetwMin *, std::string); ///< The constructor of the class.
      virtual ~Network();                    ///< The destructor of the class.

      std::list <tNode> nodes;               ///< A list of the nodes the network consists of. The nodes are described 
                                             ///< by using the sNode stucture.
      std::list <tStrand> strands;           ///< A list of all strands present in the network, described by using the 
                                             ///< sStrand structure.
      std::list <tSubCh> subchains;          ///< A list of subchains present in the network, described by using the 
                                             ///< tSubCh structure.
      std::list <tStrand *> pslip_springs;   ///< A list of slip-springs present in the network. In order to avoid
                                             ///< duplicated occurences of the slip-spring strands, we use pointers to
                                             ///< strands stored in the #nodes list defined above.
      std::vector <tBead_type> node_types;   ///< A vector which stores the desctiption of bead types present 
                                             ///< in the system.
      std::vector <tBond_type> bond_types;   ///< A vector which stores the desctiption of bond types present 
                                             ///< in the system.
      
      std::vector<std::list<tStrand *> > sorted_chains;  ///< Each element of the #sorted_chains vector consists of 
                                                         ///< a list of pointers to the internal strands a polymeric 
                                                         ///< chain consists of. The pointers refer to the array 
                                                         ///< #strands defined above.
   
      
   };
   bool pred_strand_is_end(tStrand *);    ///< A boolean function for judging whether a strand is a chain end.
}

#endif	/* NETWORK_H */

