/** @file
 *  @author Grigorios G. Megariotis (<gmegariotis@yahoo.gr>)
 *  @version 1.0
 *  @brief Random number generator based on Marsaglia's KISS (Keep it Simple and Stupid)
 *         algorithm. 
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2014 Grigorios Megariotis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef MARSAGLIA_H
#define	MARSAGLIA_H

namespace NetworkNS {

   /** @class RanMars
    *  @brief The class of the pseudorandom number generator. It is based on Marsaglia's KISS design. 
    */
   class RanMars {
   public:
      
      RanMars(int);  ///< The constructor of the random number generator class.
      
      ~RanMars();    ///< The destructor of the random number generator class.
      
      double uniform();    
      ///< Function generating unifomly distributed random numbers in [0,1).
      
      double gaussian();   
      ///< Function generating random numbers distributed according to a Gaussian distribution
      ///< centered at zero with unit standard deviation.
      
      double modified_gaussian(double mean, double stdev);
      ///< A function returning a modified Gaussian function, centered at a specified mean and 
      ///< having a pre-specified deviation.
      
      double rand_gauss(void);

      unsigned int devrand();
      
      unsigned int uint_rand();
      
   private:
      int seed;   ///< The seed used for the pseudorandom number generator.
      int save;   
      double second;
      double *u;
      int i97, j97;
      
      unsigned int x;
      unsigned int y;
      unsigned int z; 
      unsigned int c; 

   };

}


#endif	/* MARSAGLIA_H */

