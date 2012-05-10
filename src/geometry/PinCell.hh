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
class PinCell : public Object
{

public:

  typedef SP<PinCell>    SP_pincell;
  typedef Mesh::SP_mesh  SP_mesh;

  /*!
   *  \brief Constructor.
   *
   *  \param    pitch       Pin cell pitch (assumed square)
   *  \param    radii       Vector of fuel pin radii (can be zero length)
   *  \param    mat_map     Region material map (cell-center outward)
   */
  PinCell(double pitch, vec_dbl radii, vec_int mat_map);

  /// SP Constructor
  static SP<PinCell> Create(double pitch, vec_dbl radii, vec_int mat_map)
  {
    SP_pincell p;
    p = new PinCell(pitch, radii, mat_map);
    return p;
  }

  /// Return the smart pointer to my mesh.
  Mesh::SP_mesh mesh()
  {
    return d_mesh;
  }

  /// Get const reference to my mesh.
  const Mesh2D& mesh_ref() const
  {
    return *d_mesh;
  }


  /*!
   *  \brief Mesh the pin cell.
   *  \param number_meshes  Number of meshes per direction.
   *  \param flag           Switch between uniform and "best fit"
   */
  void meshify(int number_meshes, bool flag = false);

  /// DBC method.
  bool is_valid() const
  {
    return true;
  }

private:

  /// Underlying meshed object
  Mesh2D::SP_mesh d_mesh;

  /// Pin cell pitch.
  double d_pitch;

  /// Pin shell radii.
  vec_dbl d_radii;

  // Region material map.
  vec_int d_mat_map;

  /*! ======================================================================
   * @brief Determine in what region a mesh center resides.
   *
   *
   * @param  i  Horizontal index.
   * @param  j  Vertical index.
   * @return      Region index.
   */
  int find_region(int i, int j, double width);

};

} // end namespace detran

#endif /* PinCell_HH_ */

//---------------------------------------------------------------------------//
//              end of PinCell.hh
//---------------------------------------------------------------------------//
