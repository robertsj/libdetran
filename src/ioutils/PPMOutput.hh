//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PPMOutput.hh
 *  @author robertsj
 *  @date   Mar 13, 2013
 *  @brief  PPMOutput class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_ioutils_PPMOUTPUT_HH_
#define detran_ioutils_PPMOUTPUT_HH_

#include "ioutils/ioutils_export.hh"
#include "ioutils/PPMPlotter.hh"
#include "transport/State.hh"
#include "geometry/Mesh.hh"

namespace detran_ioutils
{

/**
 *  @class PPMOutput
 *  @brief Produces 2-D plots of the mesh, materials, or state.
 *
 *  This is really limited, as no color bar or other information is
 *  displayed---just colors.  It's mostly as a fast way to plot the
 *  mesh when Python+matplotlib isn't the answer.
 */
class IOUTILS_EXPORT PPMOutput
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef PPMPlotter::SP_ppmplotter         SP_ppmplotter;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_int         vec_int;
  typedef detran_geometry::Mesh::SP_mesh    SP_mesh;
  typedef detran::State::SP_state           SP_state;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //-------------------------------------------------------------------------//

  PPMOutput(SP_mesh mesh);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Initialize file
  bool initialize(const std::string filename = "detran", size_t n = 1);

  /**
   *  @brief Write a mesh map to file.
   *  @param key  Key of the map to write
   *  @return     True for successful write
   */
  bool write_mesh_map(const std::string &key);

  /**
   *  @brief Write the multigroup scalar flux moments to file.
   *  @param state  State vector container
   *  @return     True for successful write
   */
  bool write_scalar_flux(SP_state state);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Plotter
  SP_ppmplotter d_plotter;
  /// Filename base
  std::string d_filename;
  /// Mesh
  SP_mesh d_mesh;
  /// Mesh resolution
  double d_dxyx;
  /// X resolution
  size_t d_nx;
  /// Y resolution
  size_t d_ny;
  /// x values
  vec_dbl d_x;
  /// y values
  vec_dbl d_y;
  /// cells for x and y
  vec_int d_cells;

};

} // end namespace detran_ioutils

#endif /* detran_ioutils_PPMOUTPUT_HH_ */
