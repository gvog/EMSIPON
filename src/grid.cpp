/** @file
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @brief The source file containing the routines of the grid class used for the estimation of 
 *         non-bonded interactions.
 * 
 *  
 *  @section LICENSE
 *      
 *  Copyright (c) 2014 Grigorios Megariotis and Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include <cstdlib>
#include "grid.h"

namespace NetworkNS {

   Grid::Grid(double Lx, double Ly, double Lz, int nx, int ny, int nz) {

      /* Initialize the arrays with the cells: */
      ncellx = nx;
      ncelly = ny;
      ncellz = nz;
      ncells = ncellx * ncelly * ncellz;
      cells = (grid_cell*)malloc(ncells * sizeof(grid_cell));

      /*define dlx, dly, dlz and vecell: all of these are common parameters for each subcell
       *  and thus they are defined only once*/

      dlx = Lx / double(ncellx);
      idlx = 1.0 / dlx;
      dly = Ly / double(ncelly);
      idly = 1.0 / dly;
      dlz = Lz / double(ncellz);
      idlz = 1.0 / dlz;
      vcell = dlx * dly*dlz;
      ivcell = 1.0 / vcell;

      int counter = 0;
      //define Id, and the vector of each cell

      for (int k = 0; k < ncellz; k++) 

         for (int j = 0; j < ncelly; j++) 

            for (int i = 0; i < ncellx; i++) {


               cells[counter].Id = k * ncelly * ncellx + j * ncellx + i;
               cells[counter].Vec[0] = double(i + 1) * dlx;
               cells[counter].Vec[1] = double(j + 1) * dly;
               cells[counter].Vec[2] = double(k + 1) * dlz;


              
               counter = counter + 1;
            }

         

      

      //define the first neighbours of all nodes 


      int icell, jcell, kcell, klower, kupper, jlower, jupper, ilower, iupper, c_id;

      for (kcell = 0; kcell < ncellz; kcell++) {

         kupper = kcell + 1;
         if ((kcell + 1)>(ncellz - 1)) kupper = 0;
         klower = kcell - 1;
         if ((kcell - 1) < 0) klower = ncellz - 1;

         for (jcell = 0; jcell < ncelly; jcell++) {

            jupper = jcell + 1;
            if ((jcell + 1)>(ncelly - 1)) jupper = 0;
            jlower = jcell - 1;
            if ((jcell - 1) < 0) jlower = ncelly - 1;


            for (icell = 0; icell < ncellx; icell++) {

               iupper = icell + 1;
               if ((icell + 1)>(ncellx - 1)) iupper = 0;
               ilower = icell - 1;
               if ((icell - 1) < 0) ilower = ncellx - 1;

               c_id = (kcell) * ncelly * ncellx + (jcell) * ncellx + icell;

               // 1: z, y , x:
               //neigh[cell_id][0] = (kcell)*ncelly*ncellx + (jcell)*ncellx + icell;
               cells[c_id].neigh[0] = (kcell) * ncelly * ncellx + (jcell) * ncellx + icell;
               //2: same z, same y, x+
               //neigh[cell_id][1] = (kcell)*ncelly*ncellx + (jcell)*ncellx + iupper;
               cells[c_id].neigh[1] = (kcell) * ncelly * ncellx + (jcell) * ncellx + iupper;
               // 3: same z, same y, x-
               //neigh[cell_id][2] = (kcell)*ncelly*ncellx + (jcell)*ncellx + ilower;  
               cells[c_id].neigh[2] = (kcell) * ncelly * ncellx + (jcell) * ncellx + ilower;
               //4: same z, y+ , same x
               //neigh[cell_id][3] = (kcell)*ncelly*ncellx + (jupper)*ncellx + icell; 
               cells[c_id].neigh[3] = (kcell) * ncelly * ncellx + (jupper) * ncellx + icell;
               //5: same z, y- , same x
               //neigh[cell_id][4] = (kcell)*ncelly*ncellx + (jlower)*ncellx + icell;
               cells[c_id].neigh[4] = (kcell) * ncelly * ncellx + (jlower) * ncellx + icell;
               //6: same z, y+ , x+
               //neigh[cell_id][5] = (kcell)*ncelly*ncellx + (jupper)*ncellx + iupper;
               cells[c_id].neigh[5] = (kcell) * ncelly * ncellx + (jupper) * ncellx + iupper;
               //7: same z, y+ , x-
               //neigh[cell_id][6] = (kcell)*ncelly*ncellx + (jupper)*ncellx + ilower;
               cells[c_id].neigh[6] = (kcell) * ncelly * ncellx + (jupper) * ncellx + ilower;
               //8: same z, y-, x+
               //neigh[cell_id][7] = (kcell)*ncelly*ncellx + (jlower)*ncellx + iupper;
               cells[c_id].neigh[7] = (kcell) * ncelly * ncellx + (jlower) * ncellx + iupper;
               //9: same z, y-, x-
               //neigh[cell_id][8] = (kcell)*ncelly*ncellx + (jlower)*ncellx + ilower;
               cells[c_id].neigh[8] = (kcell) * ncelly * ncellx + (jlower) * ncellx + ilower;

               //10: z+1, y, x :
               //neigh[cell_id][9] = (kupper)*ncelly*ncellx + (jcell)*ncellx + icell;
               cells[c_id].neigh[9] = (kupper) * ncelly * ncellx + (jcell) * ncellx + icell;
               //11: z+1, same y, x+
               //neigh[cell_id][10] = (kupper)*ncelly*ncellx + (jcell)*ncellx + iupper;
               cells[c_id].neigh[10] = (kupper) * ncelly * ncellx + (jcell) * ncellx + iupper;
               //12: z+1, same y, x-
               //neigh[cell_id][11] = (kupper)*ncelly*ncellx + (jcell)*ncellx + ilower;
               cells[c_id].neigh[11] = (kupper) * ncelly * ncellx + (jcell) * ncellx + ilower;
               //13: z+1, y+ , same x
               //neigh[cell_id][12] = (kupper)*ncelly*ncellx + (jupper)*ncellx + icell;
               cells[c_id].neigh[12] = (kupper) * ncelly * ncellx + (jupper) * ncellx + icell;
               //14: z+1, y- , same x
               //neigh[cell_id][13] = (kupper)*ncelly*ncellx + (jlower)*ncellx + icell;
               cells[c_id].neigh[13] = (kupper) * ncelly * ncellx + (jlower) * ncellx + icell;
               //15: z+1, y+ , x+
               //neigh[cell_id][14] = (kupper)*ncelly*ncellx + (jupper)*ncellx + iupper;
               cells[c_id].neigh[14] = (kupper) * ncelly * ncellx + (jupper) * ncellx + iupper;
               //16: z+1, y+ , x-
               //neigh[cell_id][15] = (kupper)*ncelly*ncellx + (jupper)*ncellx + ilower;
               cells[c_id].neigh[15] = (kupper) * ncelly * ncellx + (jupper) * ncellx + ilower;
               //17: z+1, y-, x+
               //neigh[cell_id][16] = (kupper)*ncelly*ncellx + (jlower)*ncellx + iupper;
               cells[c_id].neigh[16] = (kupper) * ncelly * ncellx + (jlower) * ncellx + iupper;
               //18: z+1, y-, x-
               //neigh[cell_id][17] = (kupper)*ncelly*ncellx + (jlower)*ncellx + ilower;
               cells[c_id].neigh[17] = (kupper) * ncelly * ncellx + (jlower) * ncellx + ilower;

               //19: z-1, y, x :
               //neigh[cell_id][18] = (klower)*ncelly*ncellx + (jcell)*ncellx + icell;
               cells[c_id].neigh[18] = (klower) * ncelly * ncellx + (jcell) * ncellx + icell;
               //20: z-1, same y, x+
               //neigh[cell_id][19] = (klower)*ncelly*ncellx + (jcell)*ncellx + iupper;
               cells[c_id].neigh[19] = (klower) * ncelly * ncellx + (jcell) * ncellx + iupper;
               //21: z-1, same y, x-
               //neigh[cell_id][20] = (klower)*ncelly*ncellx + (jcell)*ncellx + ilower;
               cells[c_id].neigh[20] = (klower) * ncelly * ncellx + (jcell) * ncellx + ilower;
               //22: z-1, y+ , same x
               //neigh[cell_id][21] = (klower)*ncelly*ncellx + (jupper)*ncellx + icell;
               cells[c_id].neigh[21] = (klower) * ncelly * ncellx + (jupper) * ncellx + icell;
               //23: z-1, y- , same x
               //neigh[cell_id][22] = (klower)*ncelly*ncellx + (jlower)*ncellx + icell;
               cells[c_id].neigh[22] = (klower) * ncelly * ncellx + (jlower) * ncellx + icell;
               //24: z-1, y+ , x+
               //neigh[cell_id][23] = (klower)*ncelly*ncellx + (jupper)*ncellx + iupper;
               cells[c_id].neigh[23] = (klower) * ncelly * ncellx + (jupper) * ncellx + iupper;
               //25: z-1, y+ , x-
               //neigh[cell_id][24] = (klower)*ncelly*ncellx + (jupper)*ncellx + ilower;
               cells[c_id].neigh[24] = (klower) * ncelly * ncellx + (jupper) * ncellx + ilower;
               //26: z-1, y-, x+
               //neigh[cell_id][25] = (klower)*ncelly*ncellx + (jlower)*ncellx + iupper;
               cells[c_id].neigh[25] = (klower) * ncelly * ncellx + (jlower) * ncellx + iupper;
               //27: z-1, y-, x-
               //neigh[cell_id][26] = (klower)*ncelly*ncellx + (jlower)*ncellx + ilower;
               cells[c_id].neigh[26] = (klower) * ncelly * ncellx + (jlower) * ncellx + ilower;

            }

         }

      }



   }

   int Grid::find_grid_cell(const double &xnode, const double &ynode, 
                             const double &znode) {

      double xid = xnode * idlx;
      double yid = ynode * idly;
      double zid = znode * idlz;

      return(int(zid)*(ncelly)*(ncellx) + int(yid) * ncellx + int(xid));

   }
   
   Grid::~Grid(){
      free(cells);
      
      return;
   }
}