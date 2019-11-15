/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @version 3.0 (May 12, 2012)
 *  @brief  The header file containing elementary data types of bead and strand.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef _NET_TYPES_H
#define _NET_TYPES_H

#include <list>
#include <vector>

// Forward declaration, in order to create a pointer...
using namespace std;


/* Forward declaration of structure (object) strand in order to
 * create a pointer to this kind inside structure (object) node. */
struct sStrand;
struct sChain;                    //new addition, November 2012


 /**  @brief The sNode is the basic struct keeping all information relevant to a bead or a network node. 
  *          Once it is defined, it is converted to type tNode, which is used throughout the application.
  */
typedef struct sNode{
   int Id;                             ///< The identity tag of the nodal point.
   int Type;                           ///< The type of the nodal point:
                                       ///< - "1" corresponds to chain ends, 
                                       ///< - "2" corresponds to internal beads, and
                                       ///< - "3" corresponds to crosslinks.
   
   double Pos[3];                      ///< The position of the node.
   
   std:: vector <int>       OrChains;  ///< The chain or chains to which the node belongs to.
   
   std:: vector <sStrand *> pStrands;  ///< A vector of pointers to the strands the node is connected to.
   
   std:: vector <int>       SubCh;     ///< A vector of the IDs of the subchains the node is part of.
   
   std:: vector <sChain *>  pChains;   ///< A vector of pointers to the chains the node belongs to.
   
   int  node_cell;                     ///< The cell of the density estimation grid the node belongs to. 
                                       ///< This feature seems obsolete. It may be removed in future version.
   
   double n_mass;                      ///< Mass of the node in Kuhn segments
   
   double mass;                        ///< Mass of the node in g/mol (molecular weight).
   
   double r_node;                      ///< Edge length of a cube formed around the node for the estimation of 
                                       ///< nonbonded interactions.
   
   double r_star;                      ///< Edge length of a cube formed around the node for the estimation of 
                                       ///< nonbonded interactions, computed by employing a star polymer approximation.
                                       ///< This feature is obsolete. It will be removed in a future version.
}tNode;


/// @brief The sStrand is the basic abstract data type keeping all relevant information concerning a strand of a chain.
typedef struct sStrand{
   
   int Id;     ///< The identity tag of the strand.
   
   int Type;   ///< The type of the strand. 
               ///< - "1" stands for internal chain strands
               ///< - "2" stands for slip-springs
   
   int OrChain;      ///< The chain to which this strand belongs (applicable only if it is an internal strand)/
   
   bool slip_spring;          ///< A boolean variable desribing whether the strand is a slip-spring or not.
   
   unsigned int tcreation;    ///< The time when the strand was created. (It is useful for calculating its lifetime.)
   
   double spring_coeff;    ///< The spring coefficient of the strand. This quantity depends on the nature of the 
                           ///< potential used to describe the specific strand. More can be found in 
                           ///< distributions.cpp file documentation.
   
   double sq_end_to_end;   ///< The equilibrium squared end-to-end distance of the strand, i.e. 
                           ///< @f$ \left \langle R_{\rm e}^2 \right \rangle @f$.
   
   double kuhn_length;     ///< The Kuhn length of the underlying Kuhn segments of the strand, i.e. @f$ b @f$.
   
   std:: vector <tNode *> pEnds;    ///< Pointers to the nodal points the strands connects.
   
   double * pChain;                 ///< A pointer to the chain the strand belongs to. Obsolete feature.
   
} tStrand;

/// @brief A subchain of the network, treated as a vector of internal strands.
typedef std::list <tStrand> tSubCh;

/// @brief An elementary data type for reading in the information concerning a bead. 
typedef struct sBead_type{
   
   double mass;      ///< The mass of the bead type in g/mol.
   double n_mass;    ///< The mass of the bead type in number of Kuhn segments.
   double r_node;    ///< The edge length of the bead, if its mass is smeared into a cube.
} tBead_type;

/// @brief An elementary data type for reading in the information concerning a strand.
typedef struct sBond_type{
   
   double spring_coeff;    ///< The spring coefficient of the strand. This quantity depends on the nature of the 
                           ///< potential used to describe the specific strand. More can be found in 
                           ///< distributions.cpp file documentation.
   double sq_ete;          ///< The equilibrium squared end-to-end distance of the strand.
   double kuhnl;           ///< The Kuhn length of the underlying Kuhn segments of the strand, i.e. @f$ b @f$.
} tBond_type;
        
#endif	/* _NET_TYPES_H */

