/**
 * @file
 * @author  Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 * @version 3.0
 * @brief   This file contains the necessary routine for carrying out the slip-spring kinetic 
 *          MC simulation.
 * @warning This class assumes that a network has been initialized and positions of the nodes
 *          can be found at the corresponding array. 
 */

#include <iostream>
#include <list>
#include <vector>
#include <cmath>   
#include "constants.h"
#include "distributions.h"
#include "b3D_integrator.h"
#include "network.h"
#include "netmin.h"
#include "hopping.h"

#include "stdlib.h"

using namespace std;
using namespace NetworkNS;


namespace NetworkNS {

   /** The constructor takes care of opening a file to write the life-time of slip-spring to. */
   Hopping::Hopping(double hopping_rate_constant) {

      /* Open a file to write the life-time of slip-springs. */
      lifetimes_file.open("ss_lifetimes.txt", ofstream::out);
      if (not lifetimes_file.good())
         cout << "Specified slip-springs lifetime file does not exist!" << endl;

      events_file.open("events.txt", ofstream::out);
      if (not events_file.good()) 
         cout << "Specified events file does not exist!" << endl;
      
      nu_hopping_times_exp_of_barrier = hopping_rate_constant; // s^{-1} 

      return;
   }

   /** The destructor which takes care of closing the file of the slip-springs lifetimes. */
   Hopping::~Hopping(){
      /* Close the lifetimes file. */
      lifetimes_file.close();
      // gvog: Close the events file: 
      events_file.close();
      
      return;
   }
   
   
   /** @param[in] netapp A pointer to the original application.
    *  @param[in] b3D    A pointer to a Brownian Dynamics simulation scheme, in order to extract 
    *                    the current positions of the beads.
    *  @param[in] pos_array The array containing the positions of the beads.
    *  @param[in] temperature The temperature at which the rates will be calculated.
    *  @param[in] elapsed_time The time in ps elapsed from the previous call to the hopping routine.
    */
   void Hopping::hopping_step(NetwMin *netapp, const cb3D_integrator *b3D,
           double *pos_array, double temperature, double elapsed_time) {

      /** In order to develop a formalism of elementary events of slip-spring hopping, creation or 
       *  destruction, we need expressions for the rate of slippage along the chain backbone. In 
       *  order to extract the diffusivity of the slip-springs, we will proceed along the lines of 
       *  Terzis and Theodorou work. @cite Macromolecules_35_508 We describe self-diffusion along 
       *  the chain contour with the Rouse model. The Rouse model addresses the dynamics of polymers 
       *  in unentangled melts. A polymer chain is represented by a set of beads connected by 
       *  harmonic springs. The dynamics, as in our simulations, are modeled as a Brownian motion of 
       *  these tethered beads, the environment of a chain being represented as a continuum 
       *  (viscous medium), ignoring all excluded volume and hydrodynamic interactions.
       */

      /** In this model the self-diffusion of the center of the mass of the polymer is related to 
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
       */

      double dr[3], rsq;
      unsigned int transitions_performed = 0; // transitions counter
      double time_in_seconds = elapsed_time * 1.e-12; // simulation time in s
      double feng_old = 0.0;
              
#ifndef CONST_SL_SCHEME
      unsigned int new_slipsprings = 0;
#endif

      /** Hence, one must have: 
       * @f[\nu_{\rm diff} = \frac{k_{\rm B}T}{n_{\rm Kuhns/bead}b^2 \zeta} 
       * \exp{\left(-\frac{A_0}{k_{\rm B}T}\right)} @f]
       * where @f$ A_0 @f$ is a free energy per slip-spring in the equilibrium melt, which 
       * establishes a baseline for measuring free energies.
       */

      // gvog: A list of strands to be deleted at the current time step.
#ifdef CONST_SL_SCHEME
      std::list<std::pair<tStrand *, unsigned int > > to_be_deleted;
#else
      std::list<tStrand *> to_be_deleted;
#endif
      
      
      /* gvog: Loop over all slip-springs. */
      for (std::list<tStrand *>::iterator it = netapp->network->pslip_springs.begin();
              it != netapp->network->pslip_springs.end(); ++it) {

         //compute the initial free-energy of the current slip-spring
         dr[0] = pos_array[3 * ((*it)->pEnds[0]->Id - 1) + 0] 
               - pos_array[3 * ((*it)->pEnds[1]->Id - 1) + 0];
         
         dr[1] = pos_array[3 * ((*it)->pEnds[0]->Id - 1) + 1] 
               - pos_array[3 * ((*it)->pEnds[1]->Id - 1) + 1];
         
         dr[2] = pos_array[3 * ((*it)->pEnds[0]->Id - 1) + 2] 
               - pos_array[3 * ((*it)->pEnds[1]->Id - 1) + 2];
         
         // gvog: Ask for the shortest distance between the two beads:
         netapp->domain->minimum_image(dr[0], dr[1], dr[2]);

         // Calculate the free energy the spring contributes to the total free energy of the system.
#ifdef FENE_SLS         
         feng_old = e_fene    (dr, (*it)->spring_coeff, (*it)->sq_end_to_end, temperature);
#else
         feng_old = e_gaussian(dr, (*it)->spring_coeff, (*it)->sq_end_to_end, temperature);
#endif

         double cur_slipspring_sq_distance = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
         
         // gvog: Calculate the exponential of the slip-spring free energy.
         double spring_boltz_factor = exp(feng_old/boltz_const_kJoule_molK/temperature);
         // gvog: Calculate the probability in a timestep of Delta t
         double hopping_prob = nu_hopping_times_exp_of_barrier * spring_boltz_factor * time_in_seconds;
                  
         /* REMINDER:
          *  Type = 1 for chain ends.
          *  Type = 2 for internal beads.
          *  Type = 3 for crosslinks.
          */
         
         // gvog: Allocate a list of possible candidates to jump on for every end of the slip-spring:
         vector<list<pair<tNode*, double> > > sls_attchmnts(2);
         
         /* Loop over both ends of the slip-spring. */
         for (unsigned int iend = 0; iend < 2; iend++) {

            tNode* cur_sls_end = (*it)->pEnds[iend];
            
            /* Check whether the slip-spring has any of its ends connected to a chain end. */
            if (cur_sls_end->Type == 1) 
               /* The one possible point of attachment is the vacuum, so set it accordingly. */
               sls_attchmnts[iend].push_back(pair<tNode*,double>(static_cast<tNode *>(0), hopping_prob));
            
            if (cur_sls_end->Type != 3) {
               // Loop over all strands connected to the other end of the slip-spring
               for (vector<tStrand *>::iterator end_inc_strand = cur_sls_end->pStrands.begin();
                    end_inc_strand != cur_sls_end->pStrands.end(); ++end_inc_strand){
                  
                  // check which end of the incident strand we should consider:
                  if ((*end_inc_strand)->pEnds[0] == cur_sls_end) {
                     if ((*end_inc_strand)->pEnds[1]->Type != 3)
                        sls_attchmnts[iend].push_back(pair<tNode*, double>((*end_inc_strand)->pEnds[1], hopping_prob));
                  }
                  else {
                     if ((*end_inc_strand)->pEnds[0]->Type != 3)
                        sls_attchmnts[iend].push_back(pair<tNode*, double>((*end_inc_strand)->pEnds[0], hopping_prob));
                  }
               }
            }
            else /* (*it)->pEnds[0]->Type == 3 */
               cout << "#: hopping.cpp: Slip-spring attached to crosslink detected!.\n";
         }

//#define DEBUG_HOPPING
#ifdef DEBUG_HOPPING
         cout << " ---- ---- ---- ----" << endl;
         cout << "Slip - spring: " << (*it)->Id << endl;
         cout << "possible transitions for left end: " ;
         for (list<pair<tNode*, double> >::iterator isls = sls_attchmnts[0].begin(); 
              isls != sls_attchmnts[0].end(); ++isls){
            if ((*isls).first != 0)
               cout << "( " << (*isls).first->Id << ", " << (*isls).second << "), ";
            else
               cout << "( vacuum, " << (*isls).second << "), ";
         }
         cout << endl;
            
         cout << "possible transitions for right end: ";
         for (list<pair<tNode*, double> >::iterator isls = sls_attchmnts[1].begin(); 
              isls != sls_attchmnts[1].end(); ++isls){
            if ((*isls).first != 0)
               cout << "( " << (*isls).first->Id << ", " << (*isls).second << "), ";
            else
               cout << "( vacuum, " << (*isls).second << "), ";
         }
         cout << endl << " ---- ---- ---- ----" << endl << endl;
#endif

         // gvog: Check whether the sum of the transition probabilities is greater than 1.0.
         double summer = 0.0;
         for (vector<list<pair<tNode*,double > > >::iterator iend = sls_attchmnts.begin(); iend != sls_attchmnts.end(); ++iend)
            for (list<pair<tNode*, double > >::iterator end_ineigh = (*iend).begin(); end_ineigh != (*iend).end(); ++end_ineigh)
               summer += (*end_ineigh).second;

         // gvog: If the sum of the probabilities exceeded 1.0, normalize it to unity.
         if (summer > 1.0) {
            // TODO: Replace the following quick fix with something smarter:
            
            cout << "#: warning: Please change the transition probabilities. Sum = " << summer << " > 1.0\n";
            
            for (vector<list<pair<tNode*,double > > >::iterator iend = sls_attchmnts.begin(); iend != sls_attchmnts.end(); ++iend)
               for (list<pair<tNode*, double > >::iterator end_ineigh = (*iend).begin(); end_ineigh != (*iend).end(); ++end_ineigh)
                  (*end_ineigh).second /= summer;
         }

         
         bool perform_transition = false;
         double cum_prob = 0.0;
         double ran_num = netapp->my_rnd_gen->uniform();

         unsigned int end_to_jump = 0;
         tNode* target_node = 0;
         
         for (unsigned int iend = 0; iend < 2; iend++)
               for (list<pair<tNode*, double > >::iterator end_ineigh = sls_attchmnts[iend].begin(); 
                    end_ineigh != sls_attchmnts[iend].end(); ++end_ineigh)
               
               if (!perform_transition) {
                  // gvog: Add to the cumulative probability:
                  cum_prob += (*end_ineigh).second;

                  if (cum_prob > ran_num) {
                     // Here I have to store the pair (itr,jtr) and break.
                     end_to_jump = iend;
                     target_node = (*end_ineigh).first;
                     
                     perform_transition = true;
                  }
               }
         
         // gvog: In case the probability is too small, try another loop iteration.
         // gvog: pre 2016/02/02 - bug pointed by Aris Sgouros:
         // if ((cum_prob < 1.e-12) || (end_to_jump == 0 && target_node == 0))
         //    continue;
         // gvog: post 2016/02/02: no boolean condition
         
         bool destroy_slpsprng = false;

         if (perform_transition) {

            // gvog: Increase the counter of transitions performed.
            transitions_performed++;
            
            // gvog: Write the ID of the slip-spring to the file:
            events_file << (*it)->Id << " " << ((double)b3D->bd_cur_step) << endl;

#ifdef DEBUG_HOPPING
            if (target_node)
               cout << "end " << end_to_jump << " has hopped to " << target_node->Id << endl;
            else
               cout << "end " << end_to_jump << " has gone to vacuum" << endl;
            cout << "random number = " << ran_num << endl;
            
            exit(0);
#endif
            
            /* In case the attachment point exists: */
            if (target_node)
               (*it)->pEnds[end_to_jump] = target_node;
            else
               destroy_slpsprng = true;
            
            if (destroy_slpsprng)
               // gvog: The slip-spring can be deleted only if its distance is smaller than the attempt radius:
               if (cur_slipspring_sq_distance < (hopping_attempt_radius * hopping_attempt_radius))
#ifdef CONST_SL_SCHEME
                  to_be_deleted.push_back(std::pair<tStrand *, unsigned int>((*it), end_to_jump));
#else
                  to_be_deleted.push_back((*it));
#endif
         }
      }


#ifndef CONST_SL_SCHEME
      for (std::list<tStrand *>::iterator it = to_be_deleted.begin();
              it != to_be_deleted.end(); ++it) {

         // Find the element to be deleted:
         for (std::list<tStrand>::iterator jt = netapp->network->strands.begin();
                 jt != netapp->network->strands.end(); ++jt)
            if (&(*jt) == (*it)) {
               netapp->network->strands.erase(jt);
               break;
            }

         /* Write the lifetime of the slip-spring to the file: */
         lifetimes_file << (int) (b3D->bd_cur_step - (*it)->tcreation) << endl;

         // gvog: Finally, we can delete the slip-spring from the pointer array:
         netapp->network->pslip_springs.remove((*it));
      }
#else
      for (std::list<std::pair<tStrand *, unsigned int> >::iterator it = to_be_deleted.begin();
              it != to_be_deleted.end(); ++it) {

         /* 
          * gvog: In the constant number of slip-springs scheme, we have to check whether the 
          *       end-to-end distance of the spring to be destroyed is smaller than the capture
          *       radius, for the detailed balance to hold.
          */
         dr[0] = pos_array[3 * ((*it).first->pEnds[1]->Id - 1) + 0] 
               - pos_array[3 * ((*it).first->pEnds[0]->Id - 1) + 0];
         
         dr[1] = pos_array[3 * ((*it).first->pEnds[1]->Id - 1) + 1] 
               - pos_array[3 * ((*it).first->pEnds[0]->Id - 1) + 1];
         
         dr[2] = pos_array[3 * ((*it).first->pEnds[1]->Id - 1) + 2] 
               - pos_array[3 * ((*it).first->pEnds[0]->Id - 1) + 2];

         /* Apply minimum image convention to the separation vector. */
         netapp->domain->minimum_image(dr[0], dr[1], dr[2]);
         rsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

         /*
          * gvog: Only if the distance of the strand to be deleted is smaller than the tube diameter the 
          *       whole deletion/creation procedure can proceed.
          */
         if (rsq <= hopping_attempt_radius * hopping_attempt_radius) {

            // gvog: Initialize Rosenbluth weights to zero.
            double rosen_old = 0.0, rosen_new = 0.0;

            // gvog: Calculate the Rosenbluth weight in the old configuration:            
            for (std::list <tNode>::iterator jat = netapp->network->nodes.begin();
                    jat != netapp->network->nodes.end(); ++jat) {

               /*
                * gvog: We allow internal and end-end connections. Only crosslinks (Type = 3) are excluded 
                *       from possible candidates. Moreover, connections along the same chain are allowed,
                *       even with the same bead.
                */
#ifdef ALLOW_INTRAMOLECULAR
               // gvog: Here we exclude sites belonging to the same chain:
               if ((*jat).Type == 3)
#else
               if (((*jat).Type == 3) || ((*jat).OrChains[0] == (*it).first->pEnds[(*it).second]->OrChains[0]))
#endif
                  continue;

               /* check the distance from the end to the bead: */
               dr[0] = pos_array[3 * ((*jat).Id - 1) + 0] 
                     - pos_array[3 * ((*it).first->pEnds[(*it).second]->Id - 1) + 0];
               
               dr[1] = pos_array[3 * ((*jat).Id - 1) + 1] 
                     - pos_array[3 * ((*it).first->pEnds[(*it).second]->Id - 1) + 1];
                              
               dr[2] = pos_array[3 * ((*jat).Id - 1) + 2] 
                     - pos_array[3 * ((*it).first->pEnds[(*it).second]->Id - 1) + 2];
               /* Apply minimum image convention to the separation vector. */
               netapp->domain->minimum_image(dr[0], dr[1], dr[2]);
               rsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

               // gvog: Add the contribution of the strand to the Rosenbluth weight
               rosen_old += exp(-e_gaussian(dr, (*it).first->spring_coeff, (*it).first->sq_end_to_end, temperature)
                          / boltz_const_kJoule_molK / temperature);
            }

            // gvog: NEW CONFIGURATION
            // gvog: Select randomly one of the chain ends
            unsigned int nends = netapp->network->sorted_chains.size() * 2;
            unsigned int iend = (unsigned int) (netapp->my_rnd_gen->uniform() * double(nends));
            unsigned int ichain = (unsigned int) (iend / 2);

            tNode* new_end = 0;

            // gvog: Bug fix proposed by Aris Sgouros on January, 29
            if (iend % 2 == 0) 
               // gvog: Select the end of the chain:
               new_end = netapp->network->sorted_chains[ichain].back()->pEnds[1];
            else
               // gvog: Select the start of the chain:
               new_end = netapp->network->sorted_chains[ichain].front()->pEnds[0];

            /* Create a vector to store the candidates for slip-spring bridging. */
            std::vector<std::pair<tNode *, double> > candidates;
            std::pair<tNode *, double> cur_cand;

            /* Loop over all beads: */
            for (std::list <tNode>::iterator jat = netapp->network->nodes.begin();
                    jat != netapp->network->nodes.end(); ++jat) {

               // gvog: check whether the candidate is a crosslink, or belongs to the same chain.               
#ifdef ALLOW_INTRAMOLECULAR
               if ((*jat).Type == 3)
#else
               if (((*jat).Type == 3) || ((*jat).OrChains[0] == (new_end->OrChains[0])))
#endif
                  continue;

               /* check the distance from the end to the bead: */
               dr[0] = pos_array[3 * ((*jat).Id - 1) + 0] 
                     - pos_array[3 * (new_end->Id - 1) + 0];
               
               dr[1] = pos_array[3 * ((*jat).Id - 1) + 1] 
                     - pos_array[3 * (new_end->Id - 1) + 1];
               
               dr[2] = pos_array[3 * ((*jat).Id - 1) + 2] 
                     - pos_array[3 * (new_end->Id - 1) + 2];
               
               /* Apply minimum image convention to the separation vector. */
               netapp->domain->minimum_image(dr[0], dr[1], dr[2]);
               rsq = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];

               /* If the distance is smaller than the tube diameter of the polyisoprene,
                * append this neighbor to the candidates list. */
               if (rsq <= hopping_attempt_radius * hopping_attempt_radius) {
                  cur_cand.first = &(*jat);
                  cur_cand.second = e_gaussian(dr, (*it).first->spring_coeff, (*it).first->sq_end_to_end, temperature);

                  rosen_new += exp(-cur_cand.second / boltz_const_kJoule_molK / temperature);

                  candidates.push_back(cur_cand);
               }
            }

            /* Calculate the cumulative probability of each candidate: */
            double cum_prob = 0.0, ran_num = netapp->my_rnd_gen->uniform();
            tNode *sel_candidate = 0;
            for (std::vector<std::pair<tNode *, double> >::iterator icand = candidates.begin();
                    icand != candidates.end(); ++icand) {
               cum_prob += exp(-(*icand).second / boltz_const_kJoule_molK / temperature) / rosen_new;
               if (cum_prob > ran_num) {
                  sel_candidate = (*icand).first;
                  feng_new = (*icand).second;
                  break;
               }
            }

            // gvog: rosen_old -> old rosenbluth weight
            // gvog: rosen_new -> new rosenbluth weight
            // gvog: feng_old -> old energy
            // gvog: feng_new -> new energy
            double criterion =  exp((feng_new-feng_old)/boltz_const_kJoule_molK/temperature)*rosen_new/rosen_old;

            // gvog: Ask for a random number:
            ran_num = netapp->my_rnd_gen->uniform();
            if (criterion >= ran_num) {
               // gvog: The move has been accepted, update the connectivity of the slip-spring.
               (*it).first->pEnds[0] = new_end;
               (*it).first->pEnds[1] = sel_candidate;
            }

         } 
         else
            cout << "# warning: increase capture radius for the hopping scheme." << endl;
         
      }
#endif
      
      // cout << "#: hopping.cpp: slip-springs deleted = " << to_be_deleted.size() << endl;


      /** At every step of the 3D Brownian Dynamics simulation, where hopping kinetic Monte Carlo takes 
       *  place, every free end of the system can randomly create a new slip-spring with an internal 
       *  bead of a neighboring chain. This may be accomplished by a rate constant @f$k_{\rm creation} @f$.
       * The rate for the creation of a slip-spring is closely related to the probability of pairing the 
       * end ``a'' with one of its candidate mates which lie inside a sphere of prescribed radius 
       * @f$R_{\rm attempt} @f$. The definition of the probability implies that the more crowded chain ends
       * are the more probable to create a slip-spring. The number of neighbors around a chain end can be
       * tuned via the radius of the sphere within which the search takes place, @f$ R_{\rm attempt} @f$. 
       * A good estimate of @f$R_{\rm attempt} @f$ for polyisoprene (either pure or crosslinked) can be 
       * given by the tube diameter of the polymer. A computational study of the tube diameter of the 
       * polyisoprene as a function of the molecular weight has been done by the Li et al.
       * @cite Polymer_52_5867
       * The rate constant @f$k_{\rm creation}@f$ can be treated as an adjustable parameter of our model, 
       * which will be used to ensure that the average number of slip-spring present in the system is 
       * conserved throughout the simulation.
       */


#ifndef CONST_SL_SCHEME
      /* gvog:
       * Slip-spring creation event:
       * 1. loop over all chain ends of the system.
       * 2. for every chain end search in a sphere of prescribed radius for other beads
       * 3. calculate the probability of creating a new slip-spring and select one of the candidates at random
       */
      std::vector<std::list<tStrand *> >::iterator ich;
      std::list <tNode>::iterator jat;
      std::vector<tNode *> ich_ends(2);


      for (ich = netapp->network->sorted_chains.begin();
              ich != netapp->network->sorted_chains.end(); ++ich) {

         /* Extract the first and the last bead of chain ich: */
         ich_ends[0] = (*ich).front()->pEnds[0];
         ich_ends[1] = (*ich).back()->pEnds[1];


         /* Loop over both ends of the chain: */
         for (std::vector<tNode *>::iterator iend = ich_ends.begin();
                 iend != ich_ends.end(); ++iend) {

            /* Create a vector to store the candidates for slip-spring bridging. */
            std::list<tNode *> candidates;
            
            /* Loop over all beads of the system: */
            for (jat = netapp->network->nodes.begin();
                    jat != netapp->network->nodes.end(); ++jat) {

               // gvog: Do not construct a slip-spring with a crosslink:
#ifdef ALLOW_INTRAMOLECULAR
               if ((*jat).Type == 3)
#else
               if (((*jat).Type == 3) || ((*jat).OrChains[0] == (*iend)->OrChains[0]))
#endif
                  continue;

               /* check the distance from the end to the bead: */
               dr[0] = pos_array[3 * ((*jat).Id - 1) + 0] 
                     - pos_array[3 * ((*iend)->Id - 1) + 0];
               
               dr[1] = pos_array[3 * ((*jat).Id - 1) + 1] 
                     - pos_array[3 * ((*iend)->Id - 1) + 1];
               
               dr[2] = pos_array[3 * ((*jat).Id - 1) + 2] 
                     - pos_array[3 * ((*iend)->Id - 1) + 2];
               
               
               /* Apply minimum image convention to the separation vector. */
               netapp->domain->minimum_image(dr[0], dr[1], dr[2]);
               rsq = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

               /* If the distance is smaller than the tube diameter of the polyisoprene,
                * append this neighbor to the candidates list. */
               if (rsq <= hopping_attempt_radius * hopping_attempt_radius)
                  candidates.push_back(&(*jat));
            }

            // gvog: Calculate the creation probability for the current chain end:
            double creation_prob = nu_hopping_times_exp_of_barrier * (double)candidates.size() * time_in_seconds;
            
            // gvog: and check that it is lower than 1.0
            if (creation_prob > 1.0) {
               cout << "#: hopping.cpp: Please change the rate prefactor. Creation probability = " << creation_prob << " > 1.0\n";
               //exit(EXIT_FAILURE);
            }
            
            // gvog: ask for a random number to decide whether the creation should be attempted:
            double ran_num = netapp->my_rnd_gen->uniform();            
            tNode *sel_candidate = 0;
            unsigned int isel_cand = 0;
           
            if (creation_prob > ran_num){
               // gvog: we have to select one of the candidates to bridge our end with:
               ran_num = netapp->my_rnd_gen->uniform();
               isel_cand = (unsigned int)(ran_num*(double)candidates.size());
               
               unsigned int icand = 0;
               for (list<tNode *>::iterator it = candidates.begin(); it!=candidates.end(); ++it){
                  sel_candidate = (*it);
                  icand ++;
                  if (icand > isel_cand)
                     break;
               }
            }
            

            /* Create a new strand and append it the appropriate lists. */
            if (sel_candidate) {
               /* Here I have to create a new strand, and update the corresponding arrays.*/
               tStrand new_slip_spring;
               new_slip_spring.Id = netapp->network->strands.back().Id + new_slipsprings + 1;
               new_slip_spring.pEnds.resize(2);
               new_slip_spring.pEnds[0] = *iend;
               new_slip_spring.pEnds[1] = sel_candidate;
               new_slip_spring.slip_spring = true;
               new_slip_spring.OrChain = 0;
               new_slip_spring.tcreation = b3D->bd_cur_step;
               
               if (netapp->network->pslip_springs.size() > 0){
                  new_slip_spring.spring_coeff = netapp->network->pslip_springs.back()->spring_coeff;
                  new_slip_spring.sq_end_to_end = netapp->network->pslip_springs.back()->sq_end_to_end;
                  new_slip_spring.kuhn_length = netapp->network->pslip_springs.back()->kuhn_length;
               }
               else {
                  cout << "#: (warning): The first slip-spring of the system has been initialized with "
                       << "hard-coded coefficients.\n";
                  new_slip_spring.spring_coeff =  netapp->network->strands.front().spring_coeff;
                  new_slip_spring.kuhn_length = 9.58;
                  new_slip_spring.sq_end_to_end = pi_tube_diameter * pi_tube_diameter;
                  
                  /* Also, create a new bond type: */
                  tBond_type new_bond_type;
                  new_bond_type.spring_coeff = new_slip_spring.spring_coeff;
                  new_bond_type.kuhnl = new_slip_spring.kuhn_length;
                  new_bond_type.sq_ete = new_slip_spring.sq_end_to_end;
                  
                  netapp->network->bond_types.push_back(new_bond_type);
               }
               
               
               // gvog: and add it to the list of slip-springs
               netapp->network->strands.push_back(new_slip_spring);
               netapp->network->pslip_springs.push_back(&(netapp->network->strands.back()));
               new_slipsprings++;
            }
         }
      }
      
      //cout << "# hopping.cpp: new slip-springs created: " << new_slipsprings << endl;


      /* Re-count all springs in order to achieve continuous enumeration.*/
      unsigned int start_tag = netapp->network->strands.size() - netapp->network->pslip_springs.size();

      for (std::list<tStrand *>::iterator it = netapp->network->pslip_springs.begin();
              it != netapp->network->pslip_springs.end(); ++it)
         (*it)->Id = start_tag++;
#endif
      
      //cout << "#: hopping.cpp: transitions performed = " << transitions_performed - to_be_deleted.size() << endl;
      
      return;

   }

}

