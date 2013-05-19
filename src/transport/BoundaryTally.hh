//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  BoundaryTally.hh
 *  @brief BoundaryTally class definition
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_BOUNDARYTALLY_HH_
#define detran_BOUNDARYTALLY_HH_

#include "transport/transport_export.hh"
#include "transport/CoarseMesh.hh"
#include "transport/DimensionTraits.hh"
#include "transport/Equation.hh"
#include "angle/Quadrature.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran
{

/**
 *  @class BoundaryTally
 *  @brief Base class for recording functions of coarse mesh boundary fluxes
 */

template <class D>
class BoundaryTally
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum DIRECTED
  {
    X_DIRECTED,
    Y_DIRECTED,
    Z_DIRECTED
  };

  enum SENSE
  {
    NEGATIVE,
    POSITIVE
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<BoundaryTally>           SP_tally;
  typedef CoarseMesh::SP_coarsemesh                     SP_coarsemesh;
  typedef detran_angle::Quadrature::SP_quadrature       SP_quadrature;
  typedef CoarseMesh::SP_mesh                           SP_mesh;
  typedef typename EquationTraits<D>::face_flux_type    face_flux_type;
  typedef detran_utilities::size_t                      size_t;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec2_int                    vec2_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::vec3_dbl                    vec3_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param coarsemesh     Coarse mesher that has both coarse and fine meshes
   *  @param quadrature
   *  @param number_groups
   */
  BoundaryTally(SP_coarsemesh coarsemesh,
                SP_quadrature quadrature,
                const size_t  number_groups);

  /// Virtual destructor
  virtual ~BoundaryTally(){ /* ... */ }

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE --- ALL TALLIES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /**
   *  \brief Add angular flux to the current tally
   *
   *  @param  i           x mesh index
   *  @param  j           y mesh index
   *  @param  k           z mesh index
   *  @param  g           group index
   *  @param  o           octant
   *  @param  a           angle within octant
   *  @param  psi         edge angular flux
   */
  virtual void tally(const size_t i,
                     const size_t j,
                     const size_t k,
                     const size_t g,
                     const size_t o,
                     const size_t a,
                     const face_flux_type psi) = 0;

  /**
   *  @brief Add angular flux to the tally for a single incident direction
   *
   *  This is used to sweep over incident boundary containers so that
   *  the incident conditions can be avoided for the sweep.
   *
   *  @param  i           x mesh index
   *  @param  j           y mesh index
   *  @param  k           z mesh index
   *  @param  g           group index
   *  @param  o           octant
   *  @param  a           angle within octant
   *  @param  d           axis index for the incident flux
   *  @param  psi         edge angular flux
   */
  virtual void tally(const size_t i,
                     const size_t j,
                     const size_t k,
                     const size_t g,
                     const size_t o,
                     const size_t a,
                     const size_t d,
                     const double psi) = 0;

  /// Print all the partial currents (for debugging)
  virtual void display() = 0;

  /// Reset a group
  virtual void reset(const size_t group) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Coarse mesh
  SP_coarsemesh d_coarsemesh;
  /// Quadrature
  SP_quadrature d_quadrature;
  /// Number of groups
  const size_t d_number_groups;
  /// Fine mesh coarse edge flags
  vec2_int d_coarse_edge_flag;
  /// Shifts the spatial index depending on direction
  vec2_int d_octant_shift;
  /// Index of axes perpendicular to partial current
  /// \note Why does initializing require new standard?
  size_t d_perpendicular_index[3][2];

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  inline size_t index(const size_t i,
                      const size_t j,
                      const size_t k,
                      const size_t axis)
  {
    // Preconditions
    Require(axis < D::dimension);
    Require(i <= d_coarsemesh->get_coarse_mesh()->number_cells_x());
    Require(j <= d_coarsemesh->get_coarse_mesh()->number_cells_y());
    Require(k <= d_coarsemesh->get_coarse_mesh()->number_cells_z());

    size_t nx = d_coarsemesh->get_coarse_mesh()->number_cells_x();
    size_t ny = d_coarsemesh->get_coarse_mesh()->number_cells_y();
    if (axis == 0)
      return i + j * (nx + 1) + k * (nx + 1) * ny;
    else if (axis == 1)
      return i + j * nx + k * nx * (ny + 1);
    else
      return i + j * nx + k * nx * ny;

    // Postconditions ???
  }

};

//----------------------------------------------------------------------------//
template <class D>
BoundaryTally<D>::BoundaryTally(SP_coarsemesh   coarsemesh,
                                SP_quadrature   quadrature,
                                const size_t    number_groups)
  : d_coarsemesh(coarsemesh)
  , d_quadrature(quadrature)
  , d_number_groups(number_groups)
{
  // Preconditions
  Require(d_coarsemesh);
  Require(d_quadrature);

  SP_mesh mesh = d_coarsemesh->get_coarse_mesh();

  // Octant shift
  d_octant_shift.resize(3, vec_int(8, 0));
  // mu
  d_octant_shift[0][0] =  1;
  d_octant_shift[0][1] =  0;
  d_octant_shift[0][2] =  0;
  d_octant_shift[0][3] =  1;
  d_octant_shift[0][4] =  1;
  d_octant_shift[0][5] =  0;
  d_octant_shift[0][6] =  0;
  d_octant_shift[0][7] =  1;
  // eta
  d_octant_shift[1][0] =  1;
  d_octant_shift[1][1] =  1;
  d_octant_shift[1][2] =  0;
  d_octant_shift[1][3] =  0;
  d_octant_shift[1][4] =  1;
  d_octant_shift[1][5] =  1;
  d_octant_shift[1][6] =  0;
  d_octant_shift[1][7] =  0;
  // xi
  d_octant_shift[2][0] =  1;
  d_octant_shift[2][1] =  1;
  d_octant_shift[2][2] =  1;
  d_octant_shift[2][3] =  1;
  d_octant_shift[2][4] =  0;
  d_octant_shift[2][5] =  0;
  d_octant_shift[2][6] =  0;
  d_octant_shift[2][7] =  0;

  // Perpendicular indices
  d_perpendicular_index[0][0] = 1;
  d_perpendicular_index[0][1] = 2;
  d_perpendicular_index[1][0] = 0;
  d_perpendicular_index[1][1] = 2;
  d_perpendicular_index[2][0] = 0;
  d_perpendicular_index[2][1] = 1;
}

} // end namespace detran

#endif // detran_BOUNDARYTALLY_HH_

//----------------------------------------------------------------------------//
//              end of file BoundaryTally.hh
//----------------------------------------------------------------------------//
