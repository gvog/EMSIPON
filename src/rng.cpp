/**
 * @file
 * @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 * @brief  Random number generator based on Marsaglia's KISS (Keep it Simple and Stupid)
 *         algorithm. 
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Grigorios Megariotis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <cmath>
#include "rng.h"
#define PHI 0x9e3779b9


using namespace std;

namespace NetworkNS {

    unsigned int RanMars::devrand()
    {

        int fn;
        unsigned int r;
        fn = open("/dev/urandom", O_RDONLY);
        if (fn == -1)
            exit(-1); /* Failed! */
        if (read(fn, &r, 4) != 4)
            exit(-1); /* Failed! */
        close(fn);
        return r;
    }
    
    unsigned int RanMars::uint_rand()
    {
        unsigned long long t;
        x = 314527869 * x + 1234567;
        y ^= y << 5;
        y ^= y >> 7;
        y ^= y << 22;
        t = 4294584393ULL * z + c;
        c = t >> 32;
        z = t;
        return x + y + z;
    }


    
   RanMars::RanMars(int seed) {
      /* Seed variables */
        x = 123456789;
        y = 987654321;
        z = 43219876;
        c = 6543217; 
        
        x = devrand();
        while (!(y = devrand())); /* y must not be zero! */
        z = devrand();
        
        /* We don’t really need to set c as well but let's anyway... */
        /* NOTE: offset c by 1 to avoid z=c=0 */
        c = devrand() % 698769068 + 1; /* Should be less than 698769069 */

        // gvog: And warm up the generator:
        /* Also, discard the first values... */
        int i, nelements = 1000000;
        //unsigned int temp;
        for (i = 0; i < nelements; i++)
            uint_rand();

        return;
   }

   /* ---------------------------------------------------------------------- */

   RanMars::~RanMars() {
      delete [] u;
   }

   /* ----------------------------------------------------------------------
      uniform RN 
   ------------------------------------------------------------------------- */
   /// A random number generator returning number in the interval [0,1). 
   double RanMars::uniform() {
      
        double x;
        unsigned long long a;
        a = ((unsigned long long) uint_rand() << 32) + uint_rand();
        a = (a >> 12) | 0x3FF0000000000000ULL; /* Take upper 52 bits */
        *((unsigned long long *) &x) = a; /* Make a double from bits */
        
        return x - 1.0;
   }

   /* ----------------------------------------------------------------------
      gaussian RN 
   ------------------------------------------------------------------------- */

   double RanMars::gaussian() {
      double first, v1, v2, rsq, fac;

      if (!save) {
         int again = 1;
         while (again) {
            v1 = 2.0 * uniform() - 1.0;
            v2 = 2.0 * uniform() - 1.0;
            rsq = v1 * v1 + v2*v2;
            if (rsq < 1.0 && rsq != 0.0) again = 0;
         }
         fac = sqrt(-2.0 * log(rsq) / rsq);
         second = v1*fac;
         first = v2*fac;
         save = 1;
      } else {
         first = second;
         save = 0;
      }
      return first;
   }

   double RanMars::modified_gaussian(double mean, double stdev) {
      double first, v1, v2, rsq, fac;

      if (!save) {
         int again = 1;
         while (again) {
            v1 = 2.0 * uniform() - 1.0;
            v2 = 2.0 * uniform() - 1.0;
            rsq = v1 * v1 + v2*v2;
            if (rsq < 1.0 && rsq != 0.0) again = 0;
         }
         fac = sqrt(-2.0 * log(rsq) / rsq);
         second = mean + stdev * v1*fac;
         first = mean + stdev * v2*fac;
         //cout << "in_marsaglia=" << mean << endl;
         save = 1;
      } else {
         first = second;
         save = 0;
      }
      return first;
   }

   double RanMars::rand_gauss(void) {
      double v1, v2, s;

      do {
         v1 = 2.0 * uniform() - 1;
         v2 = 2.0 * uniform() - 1;

         s = v1 * v1 + v2*v2;
      } while (s >= 1.0);

      if (s == 0.0)
         return 0.0;
      else
         return (v1 * sqrt(-2.0 * log(s) / s));
   }

   
}