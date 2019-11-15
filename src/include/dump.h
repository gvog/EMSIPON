/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>) 
 *  @version 1.0 (Created on October 28, 2013)
 *  @brief The header file for trajectory file keeping.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef DUMP_H
#define	DUMP_H

#include <fstream>
#include <iostream>
#include <string>

#include "b3D_integrator.h"


namespace NetworkNS {

   /** @class dump
    *  @brief The class which dumps a snapshot of atom quantities (positions, atomic-level stresses) to one
    *         or more files every a predefined number of timesteps.
    */
   class dump {
   
   public:
   
      dump(std::string);   ///< The constructor of the dump class where the name of the dump file is specified.
      dump(const dump& orig); ///< The destructor of the class.
      void add_snapshot_to_dump(const class NetwMin *, const class cb3D_integrator *, unsigned int);
      ///< Add the current snapshot to the dump file. 
      virtual ~dump();
   
   
   private:
      std::ofstream my_file;  ///< The output file stream corresponding to the dump file.
      

   };
}
#endif	/* DUMP_H */

