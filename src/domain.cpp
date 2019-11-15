/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @brief  C++ source file containing the necessary functions for the manipulation of the 
 *          simulation domain.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2013 Georgios Vogiatzis and Grigorios Megariotis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "domain.h"

using namespace std;

namespace NetworkNS {

   Domain::Domain(std::string filename) {

      ifstream data_file(filename.c_str(), ifstream::in);
      /* Define an array of strings to hold the contents of the file. */
      std::vector <string> lines_of_file;

      for (int i = 0; i < 11; i++) {
         /* A temporary string for the current line of the file. */
         std::string current_line;
         getline(data_file, current_line);
         /* Add current line to file's array of lines. */
         lines_of_file.push_back(current_line);
      }

      /* x box dimension is stored in the 9th line: */
      string buf;
      stringstream ss;
      vector<string> tokens;

      // http://www.cplusplus.com/faq/sequences/strings/split/#boost-split

      for (int i = 0; i < 3; i++) {
         ss.flush();
         tokens.clear();
         ss << lines_of_file[8 + i];
         while (ss >> buf)
            tokens.push_back(buf);
         BoxLow[i] = atof(tokens[0].c_str());
         BoxHigh[i] = atof(tokens[1].c_str());
      }


      XBoxLen = BoxHigh[0] - BoxLow[0];
      iXBoxLen = 1.0 / XBoxLen;

      YBoxLen = BoxHigh[1] - BoxLow[1];
      iYBoxLen = 1.0 / YBoxLen;

      ZBoxLen = BoxHigh[2] - BoxLow[2];
      iZBoxLen = 1.0 / ZBoxLen;

      BoxExists = 1;
      NonPeriodic = 0;
      Xperiodic = Yperiodic = Zperiodic = 1;
      Periodicity[0] = Xperiodic;
      Periodicity[1] = Yperiodic;
      Periodicity[2] = Zperiodic;

      Boundary[0][0] = Boundary[0][1] = 0;
      Boundary[1][0] = Boundary[1][1] = 0;
      Boundary[2][0] = Boundary[2][1] = 0;

      return;
   }

   Domain::~Domain() {
      return;
   }

   void Domain::minimum_image(double &x, double &y, double &z) {
      x = x - XBoxLen * round(x * iXBoxLen);
      y = y - YBoxLen * round(y * iYBoxLen);
      z = z - ZBoxLen * round(z * iZBoxLen);
      return;
   }

   void Domain::zero_to_length_minimum_image(double &x, double&y, double &z)
 {
      x = x - XBoxLen * round(x * iXBoxLen);
      y = y - YBoxLen * round(y * iYBoxLen);
      z = z - ZBoxLen * round(z * iZBoxLen);

   }
   
   void Domain::put_in_primary_box(double *coords, int *pcoeffs){
      
      // Calculate the distance from the center of the simulation box.
      double distx = coords[0] - (BoxLow[0] + 0.5*XBoxLen);
      double disty = coords[1] - (BoxLow[1] + 0.5*YBoxLen);
      double distz = coords[2] - (BoxLow[2] + 0.5*ZBoxLen);
      
      // Calculate the periodic continuation coefficients:
      pcoeffs[0] = round(distx * iXBoxLen);
      pcoeffs[1] = round(disty * iYBoxLen);
      pcoeffs[2] = round(distz * iZBoxLen);
      
      coords[0] += (double)pcoeffs[0] * XBoxLen + BoxLow[0];
      coords[1] += (double)pcoeffs[1] * YBoxLen + BoxLow[1];
      coords[2] += (double)pcoeffs[2] * ZBoxLen + BoxLow[2];
      
      return;
   }

}