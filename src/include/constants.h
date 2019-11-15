/** 
 * @file
 * @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 * @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)  
 * @brief  Header file containing the definitions of physical constants.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Grigorios Megariotis and Georgios Vogiatzis. 
 *  All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */



#ifndef CONSTANTS_H
#define	CONSTANTS_H

double const avogadro_constant = 6.02214129e23;       ///< The Avogadro constant measured in 
                                                      ///< @f$ {\rm mol}^{-1} @f$
                                                        
double const kg_to_amus = 6.022141129e26;             ///< Conversion from kg to g/mol, 
                                                      ///< @f$ 1 {\rm kg} = @f$ 

double const amu_to_kg = 1.660538921e-27;             ///< Conversion from g/mol to kg.
                                                      ///< @f$ 1 \:{\rm g/mol}=1.660538921 
                                                      ///< \times 10^{-27} kg @f$

double const pi_monomer_mass = 68.12;                 ///< The mass of an isoprene monomer in g/mol.

double const boltz_const_Joule_molK = 8.3144621;      ///< The Boltzmann constant in J/mol/K.
double const boltz_const_kJoule_molK = 8.3144621e-3;  ///< The Boltzmann constant in kJ/mol/K
double const boltz_const_Joule_K = 1.3806488e-23;     ///< The Boltzmann constant in J/K.
double const pi_tube_diameter = 80.39;                ///< The tube diameter of polyisoprene 
                                                      ///< in @f$\text{\AA} @f$. 
                                                      ///< The estimation is based on the work of 
                                                      ///< Li et al.@cite Polymer_52_5867. 

const double c1=-5.4e02;                              ///< Constant @f$ c_1@f$ of the functional of 
                                                      ///< the nonbonded free energy, which is of 
                                                      ///< the form: @f$ \mathcal{V}_{\rm nb} 
                                                      ///< (\rho) = c_1 \rho + c_2 \rho^2 @f$.

const double c2=7.0e04;                               ///< Constant @f$ c_2@f$ of the functional 
                                                      ///< of the nonbonded free energy, which is 
                                                      ///< of the form: 
                                                      ///< @f$ \mathcal{V}_{\rm nb} (\rho) = 
                                                      ///< c_1 \rho + c_2 \rho^2 @f$.

const double tol = 1.e-4;                             ///< Tolerance for numerical comparisons.

const double PI = 3.1415926535897932384626433;        ///< The well-known @f$\mathrm{\pi} @f$ constant.

double const hopping_attempt_radius = 60.0;

#endif	/* CONSTANTS_H */

