/*
 * Assembly.hh
 *
 *  Created on: Apr 14, 2012
 *      Author: robertsj
 */

#ifndef ASSEMBLY_HH_
#define ASSEMBLY_HH_

#include "Mesh2D.hh"
#include "PinCell.hh"

namespace detran
{

class Assembly : public Mesh2D
{

public:

  typedef SP<Assembly>            SP_mesh;
  typedef Mesh2D                  Base;
  typedef PinCell::SP_mesh        SP_pincell;
  typedef std::vector<SP_pincell> vec_pincell;

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   *  \param    meshes      Number of evenly-spaced meshes per direction
   */
  Assembly(vec_pincell pincells, vec_int pincell_map);

   /*!
    *  \brief SP Constructor.
    */
   static SP<Mesh> Create(vec_pincell pincells,
                          vec_int pincell_map)
   {
     SP_mesh p;
     p = new Assembly(pincells, pincell_map);
     return p;
   }

private:

  /// Assembly pitch
  double d_pitch;

  /// Vector of SP pointers to pin cells in the assembly
  vec_pincell d_pincells;

  /// Logically 2-D map of pin cell locations
  vec_int d_pincell_map;
};

} // end namespace detran

#endif /* ASSEMBLY_HH_ */
