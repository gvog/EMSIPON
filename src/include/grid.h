/** @file
 *  @author Grigorios G. Megariotis (<gmegariotis@yahoo.gr>)
 *  @version 1.0 (July 21, 2012)
 *  @brief The class of the non-bonded energy estimation grid.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2012 Grigorios Megariotis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#ifndef GRID_H
#define	GRID_H


namespace NetworkNS {

   
   /// The sgrid_cell is the elementary struct for storing the information concerning a cell of the free energy
   /// estimation grid.
   typedef struct sgrid_cell {
      int Id;        ///< The ID of the cell.
      double Vec[3]; ///< The position of the center of the cell.
      int neigh[27]; ///< An array of the neighboring cells. The first record is the ID of the cell itself.
   } grid_cell;

   /** @class Grid
    *  @brief The nonboned free energy estimation grid class.
    */

   class Grid {
   public:

      Grid(double, double, double, int, int, int); ///< The constructor of the Grid class.
      virtual ~Grid(); ///< The destructor of the Grid class.

      double dlx; ///< The grid spacing along the @f$x@f$ direction, @f$ \Delta L_x@f$.
      double dly; ///< The grid spacing along the @f$y@f$ direction, @f$ \Delta L_y@f$.
      double dlz; ///< The grid spacing along the @f$z@f$ direction, @f$ \Delta L_z@f$.

      double idlx; 
      ///< The inverse of the grid spacing along the @f$x@f$ direction, @f$1/(\Delta L_x)@f$.
      
      double idly; 
      ///< The inverse of the grid spacing along the @f$y@f$ direction, @f$1/(\Delta L_y)@f$.
      
      double idlz; 
      ///< The inverse of the grid spacing along the @f$z@f$ direction, @f$1/(\Delta L_z)@f$.

      double vcell; ///< The volume of the grid cell, 
                    ///< @f$V_{\rm cell}=\Delta L_x\times\Delta L_y\times \Delta L_z@f$.
      double ivcell; ///< The inverse of the volume of the grid cell, i.e.@f$ 1 / V_{\rm cell} @f$.

      int find_grid_cell(const double &xnode, const double &ynode, const double &znode);
      ///< A routine for finding the cell to which a bead belongs to.

      int ncellx; ///< Number of cells along @f$x @f$ direction.
      int ncelly; ///< Number of cells along @f$y @f$ direction.
      int ncellz; ///< Number of cells along @f$z @f$ direction.
      int ncells; ///< Total number of cells.
      
      grid_cell *cells; ///< A vector containing the cells of the computation grid.
   };
}




#endif	/* GRID_H */