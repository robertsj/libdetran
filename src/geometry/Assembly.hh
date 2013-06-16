//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Assembly.hh
 *  @brief Assembly class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_ASSEMBLY_HH_
#define detran_geometry_ASSEMBLY_HH_

#include "Mesh2D.hh"
#include "PinCell.hh"

namespace detran_geometry
{

/**
 *  @class Assembly
 *  @brief Simple square assembly.
 *
 *  Assembly represents a square array of pin cells that
 *  are assumed to have the same meshing.  The meshing is
 *  also assumed to be isotropic (same in x and y) but not
 *  necessarily uniform.
 */
class GEOMETRY_EXPORT Assembly
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Assembly>    SP_assembly;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef PinCell::SP_pincell               SP_pincell;
  typedef std::vector<SP_pincell>           vec_pincell;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_utilities::vec_size_t      vec_size_t;
  typedef detran_utilities::size_t          size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    nx    Number of pins in the x-direction
   *  @param    ny    Number of pins in the y-direction
   */
  Assembly(const size_t nx, const size_t ny = 0);

  /// SP Constructor
  static SP_assembly Create(const size_t number_x,
                            const size_t number_y = 0);

  /// Return underlying meshed object.
  Mesh::SP_mesh mesh();

  /// Add a pincell
  void add_pincell(SP_pincell pin);

  /// Set the pincell map
  void set_pincell_map(const vec_int &pmap);

  vec_int pincell_map() const;

  /// Mesh the assembly.
  void finalize(vec_int pincell_map);

  /// Get dimension
  int dimension(const size_t dim = 0) const
  {
    if (dim == 0) return d_number_x;
    return d_number_y;
  }

  int number_pincells() const
  {
    return d_number_x * d_number_y;
  }

  vec_pincell pin_cells() const
  {
    return d_pincells;
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Meshed object
  Mesh2D::SP_mesh d_mesh;
  //@{
  /// Dimension per direction, e.g. 17 in 17x17.
  size_t d_number_x;
  size_t d_number_y;
  //@}
  /// Assembly pitch
  //double d_pitch;
  /// Vector of SP pointers to pin cells in the assembly
  vec_pincell d_pincells;
  /// Logically 2-D map of pin cell locations
  vec_int d_pincell_map;

};

GEOMETRY_TEMPLATE_EXPORT(detran_utilities::SP<Assembly>)
GEOMETRY_TEMPLATE_EXPORT(std::vector<detran_utilities::SP<Assembly> >)

} // end namespace detran_geometry

#endif /* detran_geometry_ASSEMBLY_HH_ */

//----------------------------------------------------------------------------//
//              end of Assembly.hh
//----------------------------------------------------------------------------//
