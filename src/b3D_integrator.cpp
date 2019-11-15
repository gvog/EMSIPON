/**
 * @file
 * @author  Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 * @author  Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 * @version 1.0
 * @brief   C++ source file implementing the Brownian Dynamics simulation of PI.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Grigorios Megariotis and Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */


#include <time.h>

#include "b3D_integrator.h"
#include "constants.h"
#include "distributions.h"
#include "hopping.h"
#include "netmin.h"



using namespace std;

namespace NetworkNS{

   cb3D_integrator::~cb3D_integrator() {
      free(bd_x);
      free(bd_x_ps);
      free(bd_mass);
      free(bd_f);
      free(bd_nb_f);
      free(bd_gamma);

      return;
   }

   /** @param[in] *init_net A pointer to the application itself, From this pointer
    *                       data concerning the network and the domain can be retrieved.
    *  @param[in] temperature The temperature at whichi the Brownian Dynamics integrator will be
    *                         initialized.
    */
   cb3D_integrator::cb3D_integrator(class NetwMin *init_net, double temperature, double slipspring_rate) {
      
      
      cur_bd_net = init_net;
      bd_cur_step = 0;
      dofs = cur_bd_net->network->nodes.size();
      dofs_3N = 3 * dofs;
      
      xshift = (double*)malloc(dofs*sizeof(double));
      yshift = (double*)malloc(dofs*sizeof(double));
      zshift = (double*)malloc(dofs*sizeof(double));
      grid_cell = (int*)malloc(dofs*sizeof(int));

      bd_x = (double*) malloc(dofs_3N * sizeof (double));
      bd_x_ps = (double*) malloc(dofs_3N * sizeof (double));
      bd_mass = (double*) malloc(dofs_3N * sizeof (double));
      bd_f = (double*) malloc(dofs_3N * sizeof (double));
      bd_f_ps = (double*) malloc(dofs_3N * sizeof (double));
      bd_nb_f = (double*) malloc(dofs_3N * sizeof (double));
      bd_gamma = (double*) malloc(dofs_3N * sizeof (double));
      
      
      bd_stress = (double**)malloc(dofs * sizeof(double*));
      for(unsigned int i = 0; i < dofs; i++)
         bd_stress[i] = (double*)malloc(6*sizeof(double));
      
      
      // Set the temperature:
      bd_temp = temperature;

      double gamma;
      unsigned int inode = 0;


      /** Klopffer et al. (\a Polymer \b 1998, \a 39, 3445 - 3449) have characterized the 
       *  rheological behavior of a series of polybutadienes and polyisoprenes over a wide
       *  range of temperatures. The viscoelastic coefficients resulting from the time-
       *  temperature superposition principle were determined. The Rouse theory modified for 
       *  undiluted polymers was used to calculate the monomeric friction coefficient, 
       *  @f$ \zeta_0@f$ from the transition zone.
       *  It was concluded that, within experimental error, a single set of WLF parameters at
       *  @f$ T_{\rm g}@f$ was adequate to characterize the relaxation dynamics irrespective of
       *  the vinyl content of the polybutadienes and polyisoprenes. 
       */

      /** The monomeric friction coefficient, @f$\zeta_0 @f$, characterizes the resistance 
       * encountered by a monomer unit moving through its surroundings. It has been shown to follow 
       * the WLF law. The variation of the monomeric friction coefficient with temperature is: 
       * @f[ \log{\zeta_0 \left(T\right)} = \log{\zeta_\infty} + \frac{C_1^{\rm g}C_2^{\rm g}}
       *    {T-T_{\rm g} + C_2^{\rm g}} @f]
       * with the parameters @f$ C_1^{\rm g} = 13.5 \pm 0.2 @f$, 
       * @f$ C_2^{\rm g} = 45 \pm 3 \;{\rm K}@f$, 
       * @f$ \log{\zeta_\infty} = -10.4 \;{\rm dyn\;s\;cm}^{-1}@f$ 
       * and @f$ T_{\rm g} = 211.15 \;{\rm K} @f$. 
       * At a temperature of 298 K, 
       * @f$\zeta_0(298\:{\rm K}) = 1.61 \times 10^{-6} \;{\rm dyn\;s\;cm}^{-1}@f$, while
       * at a temperature of  500 K,
       * @f$\zeta_0(500\:{\rm K}) = 2.63 \times 10^{-9} \;{\rm dyn\;s\;cm}^{-1}
       * = 2.63 \times 10^{-12} \;{\rm kg/s}@f$.
       */

      // Monomeric friction coefficient in  kg/s
      double monomeric_friction = 1.e-3 * pow(10.0, ((13.5 * 45.0)/(bd_temp-211.15+45)-10.4)); 
      
      /** The @f$\zeta @f$ parameter refers to a monomer of PI moving through its environment. 
       *  Since, we are dealing with larger entities, we analyze it as:
       * @f[ \zeta_0 \left(T\right) = \gamma\left(T\right) \cdot m_{\rm monomer} @f]
       * where @f$ \gamma (T)@f$ is measured in @f$ {\rm s}^{-1} @f$ and can be then multiplied by 
       * the mass of the Brownian bead which consists of several PI monomers.
       */

      // convert monomeric friction coefficient (kg/s) to g/mol/s and then divide it with the mass 
      // of a PI monomer 
      gamma = monomeric_friction / (pi_monomer_mass * amu_to_kg); // s^{-1}

      for (std::list<tNode>::iterator it = cur_bd_net->network->nodes.begin();
              it != cur_bd_net->network->nodes.end(); ++it) {

         // positions in A
         bd_x[3 * inode + 0] = (*it).Pos[0];
         bd_x[3 * inode + 1] = (*it).Pos[1];
         bd_x[3 * inode + 2] = (*it).Pos[2];

         // positions in A
         bd_x_ps[3 * inode + 0] = (*it).Pos[0];
         bd_x_ps[3 * inode + 1] = (*it).Pos[1];
         bd_x_ps[3 * inode + 2] = (*it).Pos[2];

         // mass of nodal points in g/mol
         bd_mass[3 * inode + 0] = (*it).mass;
         bd_mass[3 * inode + 1] = (*it).mass;
         bd_mass[3 * inode + 2] = (*it).mass;

         // All gamma are measured in s^{-1}
         bd_gamma[3 * inode + 0] = gamma;
         bd_gamma[3 * inode + 1] = gamma;
         bd_gamma[3 * inode + 2] = gamma;

         inode++;
      }

      
      /* Check whether all beads have the same mass and the same friction coefficient. */
      double prev_val = bd_mass[0] * bd_gamma[0];
      gamma_mass_opt = true;
      for (inode = 1; inode < dofs_3N; inode++){
         if ( ((bd_mass[inode]*bd_gamma[inode])-prev_val)
             *((bd_mass[inode]*bd_gamma[inode])-prev_val) >= tol){
            gamma_mass_opt = false;
            break;
         }
         else
            prev_val = bd_mass[inode]*bd_gamma[inode];
      }
      
      if (gamma_mass_opt) 
         cout << "#:\n#: Brownian Dynamics integrator will run in the optimized way.\n"
              << "#: --- All beads have the same mass and friction coefficient.\n" 
              << "#: --- The monomeric friction coefficient is " << monomeric_friction << " kg/s.\n"
              << "#: --- The friction coefficient is " << bd_gamma[0] << " s^{-1}.\n" 
              << "#: --- The bead friction coefficient is " << bd_mass[0]*bd_gamma[0]*amu_to_kg 
              << " kg/s.\n" 
              << "#: --- The expected Rouse diffusivity multiplied by N would be: " 
              << boltz_const_Joule_K * bd_temp / bd_mass[0] / bd_gamma[0] / amu_to_kg 
              << " m^2/s.\n#:"<< endl;
      
      // Initialize the coupled hopping scheme:
      // Here we can input the zeta parameter...
      if (cur_bd_net->network->pslip_springs.size() > 0)
         my_hopping_scheme = new Hopping(slipspring_rate);
            
      return;
   }

   
   /** @param[in] nsteps  The number of integration steps to be carried out.
    *  @param[in] dt  Integration timestep in ps.
    *  @param[in] nstout Every how many steps a report is given to the user. 
    */
   void cb3D_integrator::integrate(unsigned int nsteps, double dt, unsigned int nstout) {

      /** For large values of @f$\gamma \Delta t @f$ in the diffusive regime, when the friction is 
       *  so strong that the velocities relax within @f$\Delta t@f$. For example for 
       *  @f$\zeta_0(500\:{\rm K}) = 2.63 \times 10^{-12} \;{\rm kg/s}@f$, 
       *  @f$\gamma = 2.325 \times 10^{13}\; {\rm s}^{-1}@f$ and thus @f$\gamma \Delta t=23.25@f$ 
       *  for @f$ \Delta t = 10^{-12} \;{\rm s}@f$.
       *  The Brownian Dynamics algorithm consists of the following equation of motion: 
       *  @f[ x(t_n + \Delta t) = x(t_n) + \frac{\Delta t}{m_{bead} \gamma} 
       *  F(t_n) + X_n(\Delta t) @f]
      */  
      double current_energy = 0.0, current_nb_energy = 0.0;
   
      /** The coefficient @f$\Delta t/\left(\gamma(T) m_{\rm bead}\right)@f$ is measured in 
       * @f${\rm s^2}/{\left(\rm g/mol\right)}@f$. The forces, are measured in 
       * @f${\rm kJ}/ ({\rm mol \mathring{A}} )@f$. Thus, there is a factor of @f$ 10^{26} @f$ which
       * multiplies the force in order to make it compatible with the 
       * @f$\Delta t/\left(\gamma(T) m_{\rm bead}\right)@f$ product: 
       * @f[\frac{\rm kJ}{\rm mol \mathring{A}} = 10^{26}
       * \frac{\rm g\: \mathring{A}}{\rm mol \: s^2}@f]
       */
      density_cells = (double*) malloc(cur_bd_net->grid->ncells*sizeof(double));
      den_dx = (double*) malloc(27 * cur_bd_net->network->nodes.size() * sizeof (double));
      den_dy = (double*) malloc(27 * cur_bd_net->network->nodes.size() * sizeof (double));
      den_dz = (double*) malloc(27 * cur_bd_net->network->nodes.size() * sizeof (double));

      // Keep the time:
      tbegin = clock();
      double dt_times_inv_mass_gamma, std;
      for (unsigned int idof = 0; idof < dofs_3N; idof++)
            bd_f_ps[idof] = 0.0;

      // Ensure that there are slip-springs present in the network
      bool slipspring_hopping = false;
      if (cur_bd_net->network->pslip_springs.size() > 0) {
         slipspring_hopping = true;
         // and inform the user about the slip-spring hopping
         cout << "#: Slip-spring hopping has been enabled.\n";
      }
      else
         cout << "#: Slip-spring hopping has been disabled.\n";
      
            
      bool out_step;
      for (unsigned int istep = 0; istep < nsteps; istep++) {
         /* Update the step counter: */
         bd_cur_step ++;
         
         /* Check whether it is time to report the statistics.*/
         out_step = (istep % nstout == 0);
         
         /* Calculate bonded and non-bonded interactions.*/
         current_energy = bonded_force_calculation(out_step);
         /* Calculate non-bonded interactions every 5 timesteps.*/
         if (istep % 5 == 0)
            current_nb_energy = simpler_scheme_non_bonded_force_calculation();
                  
         /* If it's time to report, do it: */
         if (out_step)
            report(istep, current_energy, current_nb_energy);
         
         /* Keep a restart file: */
         if (istep % 20000 == 0){
            from_bd_x_to_polymer_network();
            cur_bd_net->write_network_to_lammps_data_file();
            //if (cur_bd_net->network->pslip_springs.size() > 1400)
            //   break;
         }
        
         // Hopping starts here.
         if ((istep % 1000 == 0) && slipspring_hopping)
            my_hopping_scheme->hopping_step(cur_bd_net, this, bd_x, bd_temp, 1.e3*dt);
         
                  
         /* Optimized integration scheme, in case every bead of the system has the same mass and
          * friction coefficient. */
         if (gamma_mass_opt){
            dt_times_inv_mass_gamma = (1.e-12*dt)/(bd_mass[0]*amu_to_kg*bd_gamma[0]); // s^2/kg
            std = sqrt(2.e20 * bd_temp * boltz_const_Joule_K * dt_times_inv_mass_gamma); // A
            dt_times_inv_mass_gamma *= 1.e23 / avogadro_constant; // s^2 mol / kg
            
            for (unsigned int i = 0; i < dofs_3N; i++) {  
               bd_x[i] = bd_x_ps[i] 
                       + (bd_f[i] + bd_nb_f[i]) * dt_times_inv_mass_gamma
                       + 0.5 * (bd_f[i] + bd_nb_f[i] - bd_f_ps[i]) * dt_times_inv_mass_gamma * (1.e-12*dt)
                       + cur_bd_net->my_rnd_gen->gaussian()*std;
               
               bd_x_ps[i] = bd_x[i];
               bd_f_ps[i] = bd_f[i] + bd_nb_f[i];
            }
         }
            
         else {
            
            for (unsigned int i = 0; i < dofs_3N; i++) {
               dt_times_inv_mass_gamma = (1.e-12 * dt) / (bd_mass[i]*amu_to_kg*bd_gamma[i]);//s^2/kg
               /** Random displacements, @f$X_{n} \left(\Delta t\right) @f$ are sampled from 
                * a Gaussian distribution with zero mean and width:
                * @f[ \left\langle X_n^2 \left(\Delta t\right)\right\rangle = 2 k_{\rm B}T 
                * \frac{\Delta t}{m_{\rm bead}\gamma(T)} @f]
                */
               std = sqrt(2.e20*bd_temp*boltz_const_Joule_molK*dt_times_inv_mass_gamma); // A
               dt_times_inv_mass_gamma *= 1.e23 / avogadro_constant; // s^2 mol / kg
               
               bd_x[i] = bd_x_ps[i] 
                       + (bd_f[i] + bd_nb_f[i]) * dt_times_inv_mass_gamma
                       + 0.5 * (bd_f[i] + bd_nb_f[i] - bd_f_ps[i]) * dt_times_inv_mass_gamma
                       + cur_bd_net->my_rnd_gen->gaussian()*std;
               
               bd_x_ps[i] = bd_x[i];
               bd_f_ps[i] = bd_f[i] + bd_nb_f[i];

            }
         }

      }
      

      /** Note that the full free energy of
       *  the system will include a contribution from the entropy elasticity of the strands between 
       *  nodal points, 
       */

      
      // Output the final statistics.
      current_energy = bonded_force_calculation(true);
      current_nb_energy = simpler_scheme_non_bonded_force_calculation();
      
      report(nsteps, current_energy, current_nb_energy);
      
      // Deallocate the arrays:
      free(den_dx);
      free(den_dy);
      free(den_dz);
      free(density_cells);

      // Update the network with the new positions of the beads.
      from_bd_x_to_polymer_network();
      cur_bd_net->write_network_to_lammps_data_file();

      return;
   }

   
   /** @param[in] istep The current step of the Brownian Dynamics integration.
    *  @param[in] b_energy The bonded energy of the system at the current timestep.
    *  @param[in] nb_energy The non-bonded energy of the system at the current timestep.
    */
   void cb3D_integrator::report(unsigned int istep, double b_energy, double nb_energy){

      double press_tens[6];

      // gvog: Ask for the current time:
      clock_t tend = clock();
      calculate_pressure(press_tens);

      double inv_vol = 1.0 / (  cur_bd_net->domain->XBoxLen
                              * cur_bd_net->domain->YBoxLen
                              * cur_bd_net->domain->ZBoxLen);
      
      cout << istep << "\t" << b_energy << "\t" << nb_energy << "\t"
              << cur_bd_net->network->pslip_springs.size() << "\t" 
              //<< pressure      * inv_vol << "\t" 
              << press_tens[0] * inv_vol << "\t"
              << press_tens[1] * inv_vol << "\t"
              << press_tens[2] * inv_vol << "\t"
              << press_tens[3] * inv_vol << "\t" 
              << press_tens[4] * inv_vol << "\t" 
              << press_tens[5] * inv_vol << "\t"
              << " # (" << (double) (tend - tbegin) / CLOCKS_PER_SEC << "s )" << endl;

      //cur_bd_net->my_traj_file->add_snapshot_to_dump(cur_bd_net, this, istep);

      return;
   }
   
   /// @param[out] *press_tens An array containing the six components of the box pressure tensor.
   void cb3D_integrator::calculate_pressure(double *press_tens){
      /** A function which accumulates the per-atom stresses in order to calculate the pressure
       *  of the simulation box.
       */

      // Initialize the pressure tensor to zero.
      for (unsigned int j = 0; j < 6; j++)
         press_tens[j] = 0.0;
      
      /** The per-atom stress is the negative of the per-atom pressure tensor. It is also really a 
       *  stress*volume formulation, meaning the computed quantity is in units of pressure*volume:
       * @f[\frac{\rm kJ}{\rm mol} = \frac{10^{33}}{6.022 \times 10^{23}} 
       *    \frac{\rm kg}{\rm m\: s^2} @f]
       * Thus, if the diagonal components of the per-atom stress tensor are summed for all beads in 
       * the system and the sum is divided by @f$3 V@f$, where @f$ V @f$  is the volume of the 
       * system, the result should be @f$ -p @f$, where @f$ p @f$ is the total pressure of the 
       * system.
       */ 
      
      // Accumulated the per-atom pressures to the global tensor:
      for (unsigned int i = 0; i < dofs; i++)
         for (unsigned int j = 0; j < 6; j++)
            press_tens[j] += bd_stress[i][j];
      
      // Convert per-atom pressure*vol in atm*Angstrom
      // An interesting thread concerning loop unrolling in C++: 
      // http://stackoverflow.com/questions/15275023/clang-force-loop-unroll-for-specific-loop
      for (unsigned int j = 0; j < 6; j++)
         press_tens[j] *= 1.e33 / avogadro_constant / 101.325e3;
      
      return;
   }
   
   
   
   /** A function for evaluating the forces and virials due to the bonded interactions
    *  of the system. Both entropic springs along the chain backbone and slip-springs
    *  representing entanglements are taken into account.
    *  @param[in] stress_calc Boolean variable controlling whether stress calculation will 
    *                         take place.
    */
   double cb3D_integrator::bonded_force_calculation(bool stress_calc) {

      double fenergy = 0.0; 
              
      // Initialize the forces to zero.
      for (unsigned int i = 0; i < dofs_3N; i++)
         bd_f[i] = 0.0;
      
      // Initialize the stresses to zero, if we have been asked for stress calculation:
      if (stress_calc)
         for (unsigned int i = 0; i < dofs; i++)
            for (unsigned int j = 0; j < 6; j++)
               bd_stress[i][j] = 0.0;
      

      int taga, taga3, tagb, tagb3;
      //double *sep_vec = (double*)malloc(3*sizeof(double));
      //double *grada   = (double*)malloc(3*sizeof(double));
      //double *gradb   = (double*)malloc(3*sizeof(double));

      double sep_vec[3], grada[3], gradb[3];
      
      for (std::list<tStrand>::iterator it = cur_bd_net->network->strands.begin();
              it != cur_bd_net->network->strands.end(); ++it) {
         
         /* Ask for the tags of the nodes connected to the current strand. */
         taga  = (*it).pEnds[0]->Id - 1;
         taga3 = 3*taga; 
         tagb =  (*it).pEnds[1]->Id - 1;
         tagb3 = 3*tagb;
         
         /* Form the "strand" vector, based on the x vector coming from
          * the minimizer. */
         sep_vec[0] = bd_x[tagb3 + 0] - bd_x[taga3 + 0];
         sep_vec[1] = bd_x[tagb3 + 1] - bd_x[taga3 + 1];
         sep_vec[2] = bd_x[tagb3 + 2] - bd_x[taga3 + 2];
         
         /* Apply minimum image convention. */
         cur_bd_net->domain->minimum_image(sep_vec[0], sep_vec[1], sep_vec[2]);

#ifdef FENE_SLS
         if ((*it).slip_spring)
            fenergy += f_fene( sep_vec, (*it).spring_coeff, 
                               (*it).sq_end_to_end, bd_temp, grada, gradb);
         else
#endif
         /* Calculate the spring's contribution to the free energy of the system.*/
         fenergy += f_gaussian( sep_vec, (*it).spring_coeff, 
                                (*it).sq_end_to_end, bd_temp, grada, gradb);
         
         // Accumulate the forces of the first atom:
         bd_f[taga3 + 0] += grada[0];
         bd_f[taga3 + 1] += grada[1];
         bd_f[taga3 + 2] += grada[2];
         // Accumulate the forces of the second atom:
         bd_f[tagb3 + 0] += gradb[0];
         bd_f[tagb3 + 1] += gradb[1];
         bd_f[tagb3 + 2] += gradb[2];

         /** The stress tensor of the atom @f$ i @f$, @f$\sigma_{i,\rm a,b} @f$ is given by the 
          *  following formula, where @f$\rm a@f$ and @f$\rm b @f$ take on values @f$ x @f$, 
          *  @f$ y @f$, @f$ z @f$ to generate the six components of the symmetric tensor:
          *  @f[ \sigma_{i,\rm a,b} = -\frac{1}{2}\sum_{j=1}^{N_{\rm b}(i)}
          *      \left(r_{i,a} - r_{j,a}\right)^{\rm min.im.} F_{ij,b}^{\rm min.im.}  @f]
          *  where @f$N_{\rm b}(i)@f$ stands for the number of bonds atom @f$ i @f$ participates to.
          */ 
         if (stress_calc) {
            bd_stress[taga][0] += 0.5 * sep_vec[0] * grada[0]; // xx
            bd_stress[taga][1] += 0.5 * sep_vec[1] * grada[1]; // yy
            bd_stress[taga][2] += 0.5 * sep_vec[2] * grada[2]; // zz
            bd_stress[taga][3] += 0.5 * sep_vec[0] * grada[1]; // xy
            bd_stress[taga][4] += 0.5 * sep_vec[0] * grada[2]; // xz
            bd_stress[taga][5] += 0.5 * sep_vec[1] * grada[2]; // yz

            bd_stress[tagb][0] += 0.5 * sep_vec[0] * grada[0];
            bd_stress[tagb][1] += 0.5 * sep_vec[1] * grada[1];
            bd_stress[tagb][2] += 0.5 * sep_vec[2] * grada[2];
            bd_stress[tagb][3] += 0.5 * sep_vec[0] * grada[1];
            bd_stress[tagb][4] += 0.5 * sep_vec[0] * grada[2];
            bd_stress[tagb][5] += 0.5 * sep_vec[1] * grada[2];
         }
      }

      
      
      
      return (fenergy);
   }

   
   /** The positions of the bead are updated. */
   void cb3D_integrator::from_bd_x_to_polymer_network(void) {

      unsigned int inode = 0;

      for (std::list<tNode>::iterator it = cur_bd_net->network->nodes.begin();
              it != cur_bd_net->network->nodes.end(); ++it) {
         (*it).Pos[0] = bd_x[3 * inode + 0];
         (*it).Pos[1] = bd_x[3 * inode + 1];
         (*it).Pos[2] = bd_x[3 * inode + 2];
         inode ++;
      }
      
      return;
   }

   
   
   /// A function implementing the nonbonded free energy estimation scheme. 
   double cb3D_integrator::simpler_scheme_non_bonded_force_calculation(void) {

      /** To deal with nonbonded (excluded volume and van der Waals attractive) interactions in the 
       *  network representation, we introduce a network free energy: 
       *  @f[ A_{\rm nb} = \int d^3 \mathbf{r} f\left[\rho\left(\mathbf{r}\right) \right] @f]
       * 
       *  In the above equation, @f$\rho\left(\mathbf{r}\right) @f$ is the local density (number of
       *  Kuhn segments per unit volume) at position @f$ \mathbf{r} @f$ and 
       *  @f$ f\left(\rho \right) @f$ is a free energy density (free energy per unit volume). 
       *  Expressions for @f$f\left(\rho \right) @f$ may be extracted from an equation of state.
       *  Here the plan is to invoke a simple expression for @f$f\left(\rho \right) @f$, in the 
       *  form of a Taylor expansion, 
       *  @f[ f\left(\rho \right) = C\rho + B\rho^2 @f]
       *  with @f$ C @f$, @f$ B @f$ fitted such that the volumetric properties (pressure and 
       *  compressibility at mean density of interest) are reproduced.   
       */ 
       
      /** Local density will be resolved only at the level of entire cells, defined by passing an 
       *  orthogonal grid through the entire system. The free energy of the system is approximated 
       *  by@f[ A_{\rm nb} = \sum_{\rm cells} V_{\rm cell}^{\rm acc} 
       *  f\left( \rho_{\rm cell}\right) @f]
       *  where @f$ V_{\rm cell}^{\rm acc} @f$, the accessible of a cell, is the volume of the 
       *  rectangular parallelepiped defining the cell minus the volume of any parts of 
       *  nanoparticles that may find themselves in the cell. 
       */
      
      /** The cell density @f$\rho_{\rm cell} @f$ must be defined based on the nodal
       *  points in and around the cell, each nodal point contributing a mass equal to the node's 
       *  mass. Each nodal point @f$ j @f$ has mass @f$ n_j @f$ (in Kuhn segments) and a 
       *  characteristic size @f$ R_j @f$. We will discuss below how these quantities depend on the 
       *  node's molecular characteristics. We denote the position vector of node @f$ j @f$ by 
       *  @f$ \mathbf{r}_j = \left( x_j, y_j, z_j \right) @f$.
       *  The cell dimensions along the @f$ x @f$, @f$ y @f$, @f$ z @f$ directions will be denoted 
       *  as @f$ L_x @f$, @f$ L_y @f$, @f$ L_z @f$, respectively.
       */
      
      /** We will focus on a cell extending between @f$ x_{\rm cell}-L_x @f$ and @f$ x_{\rm cell} @f$
       *  along the @f$ x @f$-direction, between @f$ y_{\rm cell} - L_y@f$ and @f$ y_{\rm cell} @f$
       *  along the @f$ y @f$-direction, and between @f$ z_{\rm cell} - L_z @f$ and @f$ z_{\rm cell} @f$
       *  along the @f$ z @f$-direction. In the regular grid considered, if @f$(0,0,0)@f$ is taken as 
       *  one of the grid points @f$ x_{\rm cell} @f$, @f$ y_{\rm cell} @f$, and @f$ z_{\rm cell} @f$
       *  will be integer multiples of @f$ L_x @f$, @f$ L_y @f$ and @f$ L_z @f$, respectively.
       */
      
      /** In the following we will assume that 
       *  @f[ R_j < \min{ \left(L_x, L_y, L_z \right) } @f]
       */
      
      /** The simplest option for relating the positions and masses of the node to 
       *  @f$ \rho_{\rm cell} @f$ is to envision each node @f$ j @f$ as a cube containing @f$ n_j @f$
       *  Kuhn segments, of edge length @f$ R_j @f$, centered at @f$ \mathbf{r}_j @f$.
       *  Node @f$ j @f$ will contribute to a cell if its cube (cube @f$ j @f$) overlaps with the cell. 
       *  Note that, for this to happen, it is not necessary that the nodal position of the center, 
       *  @f$ \mathbf{r}_j @f$, lie in the cell. 
       *  The mass (number of Kuhn segments) contributed by the node to the cell is:
       *  @f[ n_{j,{\rm cell}} = n_j \frac{V_{{\rm cube}\; j\cap{\rm cell}}}{V_{{\rm cube}\;j}} @f]
       *  with @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ being the volume of the intersection of cube 
       *  @f$ j @f$, associated with node @f$ j @f$, and the considered cell, while 
       *  @f$ V_{{\rm cube}\;j} = R_j^3 @f$ is the volume of cube @f$ j @f$.
       */ 
      
 
      // Variables holding the volume of a cell.
      double vcube_cell, vx, vy, vz; 
      
      // The nonbonded contribution to the free energy of the network.
      double f_nb_energy = 0.0;

      //at this point we may need to call a subroutine that will convert bd_x elements to positions

      int l; // indices used to find the parent cell and its first neighbours
      
      for (int i = 0; i < cur_bd_net->grid->ncells; i++)
         density_cells[i] = 0.0;

      int cur_node = 0, i, j, cur_elem, max_node = cur_bd_net->network->nodes.size();
      double dx, dy, dz, xl, yl, zl, half_rnode, mass_over_rnode3;

      half_rnode = 0.5 * cur_bd_net->network->nodes.front().r_node;
      mass_over_rnode3 = cur_bd_net->network->nodes.front().n_mass
              / (cur_bd_net->network->nodes.front().r_node
              * cur_bd_net->network->nodes.front().r_node
              * cur_bd_net->network->nodes.front().r_node);
      

      for (cur_node = 0; cur_node < max_node; cur_node++) {
         /*loop over the (*it).node_cell itself and its first neighbors (it is always equal 
          * to 27, or 26 starting the numbering from zero*/
         /*expressed in Angstrom^3. node coordinates have to be shifted
          * boxl/2.0 so as to be embedded into a grid extended from
          * zero to cur_bd_net->domain->XBoxLen */

         //half_rnode = 0.5 * (*it).r_node;
         //mass_over_rnode3 = (*it).n_mass / (*it).r_node / (*it).r_node / (*it).r_node;

         xshift[cur_node] = bd_x[3 * cur_node];
         yshift[cur_node] = bd_x[3 * cur_node + 1];
         zshift[cur_node] = bd_x[3 * cur_node + 2];

         //return the nodes back into the primary box 
         cur_bd_net->domain->minimum_image(xshift[cur_node], yshift[cur_node], zshift[cur_node]);

         //shift the node position to a shifted simulation box that contains the  grid   
         xshift[cur_node] += 0.5 * cur_bd_net->domain->XBoxLen;
         yshift[cur_node] += 0.5 * cur_bd_net->domain->YBoxLen;
         zshift[cur_node] += 0.5 * cur_bd_net->domain->ZBoxLen;

         grid_cell[cur_node] = cur_bd_net->grid->find_grid_cell(xshift[cur_node], yshift[cur_node], 
                                          zshift[cur_node]);


         for (j = 0; j < 27; j++) {
            // find the neighbours of (*it).node_cell, zero corresponds to the cell itself
   
            l = cur_bd_net->grid->cells[grid_cell[cur_node]].neigh[j];

            /*find the intersection of cube formed by node (*it).Id, how to define 
             * cell_vec[l][0:2]: vector of cell l, is used for the calculation of 
             * vcube_cell,*/

            // if statements for vx
            // for the computation of minimum images of xshift, yshift, zshift with respect 
            // to xcell, ycell and zcell respectively
            dx = xshift[cur_node] - cur_bd_net->grid->cells[l].Vec[0];
            dy = yshift[cur_node] - cur_bd_net->grid->cells[l].Vec[1];
            dz = zshift[cur_node] - cur_bd_net->grid->cells[l].Vec[2];

            // minimum images in a box from zero to box_l
            //where xl, yl and zl updated values, due to minimum image, of xshift, yshift and zshift
            cur_bd_net->domain->zero_to_length_minimum_image(dx, dy, dz);

            xl = cur_bd_net->grid->cells[l].Vec[0] + dx;
            yl = cur_bd_net->grid->cells[l].Vec[1] + dy;
            zl = cur_bd_net->grid->cells[l].Vec[2] + dz;


            /** Under the condition @f$ R_j < \min{ \left(L_x, L_y, L_z \right) } @f$, 
             *  @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ is obtainable as:
             * \f{eqnarray*}{
             *  V_{{\rm cube}\; j\cap{\rm cell}} & = & 
             *    \max{\left\{\left[ \min{\left(x_j + \frac{R_j}{2}, x_{\rm cell} \right)} 
             *                 -\max{\left(x_j - \frac{R_j}{2}, x_{\rm cell} -L_x \right)}\right] 
             *         , 0 \right\}} \\
             * & \times & 
             *    \max{\left\{\left[ \min{\left(y_j + \frac{R_j}{2}, y_{\rm cell} \right)} 
             *                  -\max{\left(y_j - \frac{R_j}{2}, y_{\rm cell} -L_y \right)}\right] 
             *         , 0 \right\}} \\
             * & \times & 
             *    \max{\left\{\left[ \min{\left(z_j + \frac{R_j}{2}, z_{\rm cell} \right)} 
             *                  -\max{\left(z_j - \frac{R_j}{2}, z_{\rm cell} -L_z \right)}\right] 
             *         , 0 \right\}} 
             * \f}
             */
            vx = max(min(xl + half_rnode, cur_bd_net->grid->cells[l].Vec[0])
               - max(xl-half_rnode, cur_bd_net->grid->cells[l].Vec[0]-cur_bd_net->grid->dlx), 0.0);
            
            vy = max(min(yl + half_rnode, cur_bd_net->grid->cells[l].Vec[1])
               - max(yl-half_rnode, cur_bd_net->grid->cells[l].Vec[1]-cur_bd_net->grid->dly), 0.0);
            
            vz = max(min(zl + half_rnode, cur_bd_net->grid->cells[l].Vec[2])
               - max(zl-half_rnode, cur_bd_net->grid->cells[l].Vec[2]-cur_bd_net->grid->dlz), 0.0);

            vcube_cell = vx * vy * vz;

            /** As defined by the above equation, @f$ V_{{\rm cube}\; j\cap{\rm cell}} @f$ is a 
             *  linear function of the node coordinates. Clearly, if cube @f$ j @f$ lies entirely 
             *  within the cell, @f$ V_{{\rm cube}\; j\cap{\rm cell}} = V_{{\rm cube}\:j} @f$ and, 
             *  consequently, @f$ n_{j, {\rm cell}} = n_j @f$. If however, the borders of cube 
             *  @f$ j @f$ intersect the borders of the considered cell, then node @f$ j @f$ will 
             *  contribute a mass @f$ n_{j,{\rm cell}} < n_{j} @f$ to the cell. The total mass 
             *  contributed by bead @f$ j @f$ to all cells in which it participates will always 
             *  be @f$n_{j}@f$.
             */
            
            cur_elem = 27 * cur_node + j;

            if ((xl > cur_bd_net->grid->cells[l].Vec[0] - cur_bd_net->grid->dlx - half_rnode) &&
                    (xl < cur_bd_net->grid->cells[l].Vec[0] - cur_bd_net->grid->dlx + half_rnode))
               den_dx[cur_elem] = mass_over_rnode3 * vy * vz / cur_bd_net->grid->vcell;
            else if ((xl > cur_bd_net->grid->cells[l].Vec[0] - half_rnode) &&
                    (xl < cur_bd_net->grid->cells[l].Vec[0] + half_rnode))
               den_dx[cur_elem] = -mass_over_rnode3 * vy * vz / cur_bd_net->grid->vcell;
            else
               den_dx[cur_elem] = 0.0;

            if ((yl > cur_bd_net->grid->cells[l].Vec[1] - cur_bd_net->grid->dly - half_rnode) &&
                    (yl < cur_bd_net->grid->cells[l].Vec[1] - cur_bd_net->grid->dly + half_rnode))
               den_dy[cur_elem] = mass_over_rnode3 * vx * vz / cur_bd_net->grid->vcell;
            else if ((yl > cur_bd_net->grid->cells[l].Vec[1] - half_rnode) &&
                    (yl < cur_bd_net->grid->cells[l].Vec[1] + half_rnode))
               den_dy[cur_elem] = -mass_over_rnode3 * vx * vz / cur_bd_net->grid->vcell;
            else
               den_dy[cur_elem] = 0.0;

            if ((zl > cur_bd_net->grid->cells[l].Vec[2] - cur_bd_net->grid->dlz - half_rnode) &&
                    (zl < cur_bd_net->grid->cells[l].Vec[2] - cur_bd_net->grid->dlz + half_rnode))
               den_dz[cur_elem] = mass_over_rnode3 * vx * vy / cur_bd_net->grid->vcell;
            else if ((zl > cur_bd_net->grid->cells[l].Vec[2] - half_rnode) &&
                    (zl < cur_bd_net->grid->cells[l].Vec[2] + half_rnode))
               den_dz[cur_elem] = -mass_over_rnode3 * vx * vy / cur_bd_net->grid->vcell;
            else
               den_dz[cur_elem] = 0.0;

            /** The density @f$ \rho_{\rm cell} @f$ in the considered cell is estimated as: 
             *  @f[ \rho_{\rm cell} = \frac{1}{V_{\rm cell}^{\rm acc}} \sum_j n_{j,{\rm cell}} @f]
             *  Clearly, only nodal points @f$ j @f$ whose cubes have a nonzero overlap with the 
             *  considered cell will contribute to the above summation. The positions vectors 
             *  @f$ \mathbf{r}_j @f$ of these beads will necessarily lie within the considered cell
             *  or its immediate neighbors. 
             */
            
            density_cells[l] += mass_over_rnode3*vcube_cell;
         }
      }

      /** 
       * The non-bonded energy is considered to be a quadratic function of the density, 
       * i.e., 
       * @f[ A_{\rm nb} = \sum_{i\in \:{\rm cells}} V_{{\rm cell},i} 
       * \left ( C_1\rho_i +C_2 \rho_i^2 \right)
       * @f]
       */
      
      for (i = 0; i < cur_bd_net->grid->ncells; i++) {
         density_cells[i] *= cur_bd_net->grid->ivcell; //expressed in kuhn segments per Angstom^3 
         f_nb_energy += cur_bd_net->grid->vcell 
                      * (c1 * density_cells[i] + c2 * density_cells[i] * density_cells[i]);
      }

      
      /** The precise conditions for cube @f$ j @f$ to have common points with the considered 
       *  cell are:
       *  @f[ x_{\rm cell} - L_x < x_j + \frac{R_j}{2} < x_{\rm cell} + R_j @f]
       *  @f[ y_{\rm cell} - L_y < y_j + \frac{R_j}{2} < y_{\rm cell} + R_j @f]
       *  @f[ z_{\rm cell} - L_z < z_j + \frac{R_j}{2} < z_{\rm cell} + R_j @f]
       */
      

      // new nested for-loops for updating the 3N vector of derivatives 
      double fx, fy, fz; // force components on a nod due to non-bonded interactions
      // used for updating the array of derivatives[3N] 
      for (cur_node = 0; cur_node < max_node; cur_node++){
         fx = 0.0;
         fy = 0.0;
         fz = 0.0;
      
         /** According to the above approach, the force on node @f$ j @f$ due to nonbonded 
          *  interactions is:
          *  @f[ \mathbf{F}_j = - \nabla_{\mathbf{r}_j} A_{\rm nb} 
          *  = - \sum_{\substack{{\rm cells\; having\; common} \\ {\rm points\: with \: cube\;} j}}
          *  V_{\rm cell}^{\rm acc} \left. \frac{d f}{d \rho}\right|_{\rho = \rho{\rm cell}} 
          *  \nabla_{\mathbf{r}_j} \rho_{\rm cell}
          *  @f]
          */
         
         for (j = 0; j < 27; j++) {

            cur_elem = 27 * cur_node + j;
         
            l = cur_bd_net->grid->cells[grid_cell[cur_node]].neigh[j];

            fx -= cur_bd_net->grid->vcell * (c1 + 2.0 * c2 * density_cells[l]) * den_dx[cur_elem];
            fy -= cur_bd_net->grid->vcell * (c1 + 2.0 * c2 * density_cells[l]) * den_dy[cur_elem];
            fz -= cur_bd_net->grid->vcell * (c1 + 2.0 * c2 * density_cells[l]) * den_dz[cur_elem];
         }

         bd_nb_f[3 * cur_node]     = fx;
         bd_nb_f[3 * cur_node + 1] = fy;
         bd_nb_f[3 * cur_node + 2] = fz;

      }

      return (f_nb_energy);
   }
   
   
   
   void cb3D_integrator::cell_density_nodal_points() {

      double xnew, ynew, znew;
      int Id;
      int *hist_grid = (int*) malloc(cur_bd_net->grid->ncells * sizeof (int));

      for (int i = 0; i < cur_bd_net->grid->ncells; i++)
         hist_grid[i] = 0;

      for (unsigned int inode = 0; inode < dofs; inode++) {

         xnew = bd_x[3 * inode + 0];
         ynew = bd_x[3 * inode + 1];
         znew = bd_x[3 * inode + 2];
         cur_bd_net->domain->minimum_image(xnew, ynew, znew);
         xnew = xnew + cur_bd_net->domain->XBoxLen / 2.0;
         ynew = ynew + cur_bd_net->domain->YBoxLen / 2.0;
         znew = znew + cur_bd_net->domain->ZBoxLen / 2.0;
         Id = cur_bd_net->grid->find_grid_cell(xnew, ynew, znew);
         hist_grid[Id] = hist_grid[Id] + 1;

      }

      for (int i = 0; i < cur_bd_net->grid->ncells; i++)
         fprintf(p_cell_density, "%d " " %d\n", i, hist_grid[i]);

      return;
   }

}
