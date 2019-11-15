/** @file
 *  @author Grigorios Megariotis
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @brief  Header file accompanying the ``distributions.cpp'' C++ source file.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2014 Grigorios Megariotis and Georgios Vogiatzis. 
 *  All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef _DISTRIBUTIONS_H_
#define _DISTRIBUTIONS_H_

double f_gaussian(const double *, const double &, const double &, const double &, double *, double *);
double e_gaussian(const double *, const double &, const double &, const double &);

#endif
