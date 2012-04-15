//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PinCell.hh
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PinCell_HH_
#define PinCell_HH_

// Detran
#include "Mesh2D.hh"

// Utilities
//#include "Definitions.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class PinCell
 *  \brief Two-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class PinCell : public Mesh2D
{

public:

  typedef SP<PinCell>   SP_mesh;
  typedef Mesh2D        Base;

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   *  \param    meshes      Number of evenly-spaced meshes per direction
   */
  PinCell(double pitch, vec_dbl radii, vec_int mat_map, int meshes);

   /*!
    *  \brief SP Constructor.
    */
   static SP<Mesh> Create(double pitch,
                          vec_dbl radii,
                          vec_int mat_map,
                          int meshes)
   {
     SP_mesh p;
     p = new PinCell(pitch, radii, mat_map, meshes);
     return p;
   }

private:

    double d_pitch;

    vec_dbl d_radii;

    /*! ======================================================================
     * @brief Determine in what region a mesh center resides.
     *
     *
     * @param  i  Horizontal index.
     * @param  j  Vertical index.
     * @return      Region index.
     */
    int find_region(int i, int j);

};

} // end namespace detran

#endif /* PinCell_HH_ */

//---------------------------------------------------------------------------//
//              end of PinCell.hh
//---------------------------------------------------------------------------//
