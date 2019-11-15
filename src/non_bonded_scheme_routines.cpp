/** @file
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @brief The C++ source file containing the functions of the smoothed non-bonded free 
 *         energy estimation scheme.
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

void phi_function(double x, double xj, double Rj, double delta, double &phi_value) {

   double inv_delta = 1.0 / delta;

   if (x <= xj - 0.5 * Rj - 0.5 * delta)
      phi_value = 0.0;

   else if (x > xj - 0.5 * Rj - 0.5 * delta && x <= xj - 0.5 * Rj + 0.5 * delta) {
      double temp = 2.0 * (x - xj) * inv_delta + Rj*inv_delta;
      double tempsq = temp*temp;

      phi_value = (8.0 + 15.0 * temp 
                - 10.0 * tempsq * temp + 3.0 * tempsq * tempsq * temp) / 16.0 / Rj;
   } else if (x > xj - 0.5 * Rj + 0.5 * delta && x <= xj + 0.5 * Rj - 0.5 * delta)
      phi_value = 1.0 / Rj;

   else if (x > xj + 0.5 * Rj - 0.5 * delta && x <= xj + 0.5 * Rj + 0.5 * delta) {
      double temp = 2.0 * (xj - x) * inv_delta + Rj*inv_delta;
      double tempsq = temp*temp;

      phi_value = (8.0 + 15.0 * temp 
                  - 10.0 * tempsq * temp + 3.0 * tempsq * tempsq * temp) / 16.0 / Rj;
   }
   else if (x > xj + 0.5 * Rj + 0.5 * delta)
      phi_value = 0.0;
}

void integral_minus_infinity_to_x(double xasterisk, double xj, double Rj, double delta,
        double &integral_value) {

   if (xasterisk <= xj - Rj / 2.0 - delta / 2.0)
      integral_value = 0.0;

   else if (xasterisk > xj - Rj / 2.0 - delta / 2.0 && xasterisk <= xj - Rj / 2.0 + delta / 2.0)
      integral_value = delta * (5.0 / 4.0 + 4.0 * (2.0 * (xasterisk - xj) / delta + Rj / delta)
           + 15.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 2)\
                  - 5.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 4)
           + 1.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 6)) / 16.0 / Rj;

   else if (xasterisk > xj - Rj / 2.0 + delta / 2.0 && xasterisk <= xj + Rj / 2.0 - delta / 2.0)
      integral_value = (xasterisk - xj) / Rj + 1.0 / 2.0;

   else if (xasterisk > xj + Rj / 2.0 - delta / 2.0 && xasterisk <= xj + Rj / 2.0 + delta / 2.0)
      integral_value = 1.0 - delta * (5.0 / 4.0 + 4.0 * (2.0 * (xj - xasterisk) / delta + Rj / delta) 
                     + 15.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 2)
                     - 5.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 4) 
                     + 1.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 6)) / 16.0 / Rj;

   else if (xasterisk > xj + Rj / 2.0 + delta / 2.0)
      integral_value = 1.0;

}

void integral_x_to_plus_infinity(double xasterisk, double xj, double Rj, 
                                 double delta, double &integral_value) {

   if (xasterisk <= xj - Rj / 2.0 - delta / 2.0)
      integral_value = 1.0;

   else if (xasterisk > xj - Rj / 2.0 - delta / 2.0 && xasterisk <= xj - Rj / 2.0 + delta / 2.0)
      integral_value = 1.0 - delta * (5.0 / 4.0 + 4.0 * (2.0 * (xasterisk - xj) / delta + Rj / delta) 
                     + 15.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 2)
                     - 5.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 4) 
                     + 1.0 / 4.0 * pow(2.0 * (xasterisk - xj) / delta + Rj / delta, 6)) / 16.0 / Rj;

   else if (xasterisk > xj - Rj / 2.0 + delta / 2.0 && xasterisk <= xj + Rj / 2.0 - delta / 2.0)
      integral_value = (xj - xasterisk) / Rj + 1.0 / 2.0;

   else if (xasterisk > xj + Rj / 2.0 - delta / 2.0 && xasterisk <= xj + Rj / 2.0 + delta / 2.0)
      integral_value = delta * (5.0 / 4.0 + 4.0 * (2.0 * (xj - xasterisk) / delta + Rj / delta) 
                     + 15.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 2)
                     - 5.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 4) 
                     + 1.0 / 4.0 * pow(2.0 * (xj - xasterisk) / delta + Rj / delta, 6)) / 16.0 / Rj;

   else if (xasterisk > xj + Rj / 2.0 + delta / 2.0)
      integral_value = 0.0;

}
