//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER2D_HH_
#define SWEEPER2D_HH_

// Other libtran headers
#include "SweeperBase.hh"
#include "Equation_DD_2D.hh"

// libtran utilities
#include "Definitions.hh"
#include "GenException.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Sweeper2D
 * \brief Sweeper for 2-D discrete ordinates problems.
 *
 * For a 2D mesh, the sweep must go along one axis first.  Here, we
 * sweep along the x axis first.  The mesh cells are considered in terms of
 * their horizontal and vertical edges.  An update vertical edge angular
 * flux is propagated immediately to its neighbor, while the horizontal
 * edge flux must be saved for the next row of cells.
 *
 */
//---------------------------------------------------------------------------//
class Sweeper2D : public SweeperBase
{

public:

  typedef detran_utils::SP<Sweeper2D>       SP_sweeper;
  typedef State::SP_state                   SP_state;
  typedef detran_utils::InputDB::SP_input   SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;
  typedef detran_utils::vec_dbl             edge_flux_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   */
  Sweeper2D(SP_input input,
            SP_mesh mesh,
            SP_material material,
            SP_quadrature quadrature,
            SP_state state);

  /*!
   *  \brief Sweep over all angles and space.
   *
   */
  inline void sweep(moments_type &phi,
                    moments_type &source);

private:

  /// Horizontal face fluxes
  edge_flux_type d_psi_h;

  /// Vertical face fluxes
  edge_flux_type d_psi_v;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Sweeper2D.i.hh"

#endif /* SWEEPER2D_HH_ */
