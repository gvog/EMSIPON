/** 
 *  @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @version 1.0 (May 16, 2012)
 *  @brief Header file containing the definitions of Domain class.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef DOMAIN_H
#define	DOMAIN_H

#include <string>

namespace NetworkNS {
   
   /** @class Domain
    *  @brief The class of the simulation domain.
    */
   
   class Domain{
      public:
         int BoxExists;       ///< An integer variable denoting whether the simulation box exists or no.
                              ///< 0 = not yet created, 1 = exists
         int NonPeriodic;     ///< An integer variable denoting whether the simulation box is periodic or no.
                              ///< 0 = periodic in all 3 dims
                              ///< 1 = periodic or fixed in all 6
                              ///< 2 = shrink-wrap in any of 6
         
         int Xperiodic;       ///< Periodicity along @f$ x @f$ direction, 0 = non-periodic, 1 = periodic.
         int Yperiodic;       ///< Periodicity along @f$ y @f$ direction, 0 = non-periodic, 1 = periodic.
         int Zperiodic;       ///< Periodicity along @f$ z @f$ direction, 0 = non-periodic, 1 = periodic.
         
         int Periodicity[3];  ///< xyz periodicity as an array.
         
         int Boundary[3][2];  ///< Settings for 6 boundaries
                              ///< 0 = periodic, 1 = fixed non-periodic, 2 = shrink-wrap non-periodic
                              ///< 3 = shrink-wrap.
         
         double BoxLow[3];    ///< Orthogonal box global lower bounds along all three directions.
         double BoxHigh[3];   ///< Orthogonal box global lower bounds along all three directions.
         
         
         double XBoxLen;      ///< Simulation box edge length along @f$ x @f$ direction, @f$L_x@f$.
         double YBoxLen;      ///< Simulation box edge length along @f$ y @f$ direction, @f$L_y@f$.
         double ZBoxLen;      ///< Simulation box edge length along @f$ z @f$ direction, @f$L_z@f$.
         
         double iXBoxLen;     ///< Inverse simulation box edge length along @f$ x @f$ direction, @f$ 1/ L_x@f$.
         double iYBoxLen;     ///< Inverse simulation box edge length along @f$ y @f$ direction, @f$ 1/ L_y@f$.
         double iZBoxLen;     ///< Inverse simulation box edge length along @f$ z @f$ direction, @f$ 1/ L_z@f$.
         
         
         Domain(std::string); ///< The constructor of the Domain class, where only the name of the data file 
                              ///< has to be provided.
         
         virtual ~Domain();   ///< The destructor of the Domain class.
         
         void minimum_image(double &x, double &y, double &z);  
         ///< Apply minimum image convention along the three spatial directions.
         
         void zero_to_length_minimum_image(double &, double&, double &);
         ///< Put the coordinates inside the primary simulation box.
         
         void put_in_primary_box(double *, int *);
         ///< Put the coordinates inside the primary simulation box.

   };
}

#endif	/* DOMAIN_H */

