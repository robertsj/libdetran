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

// Libtran
#include "SweeperBase.hh"
#include "SweepSource.hh"
#include "Equation_DD_2D.hh"

// Utilities
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
class Sweeper2D : public SweeperBase<_2D>
{

public:

  enum FACE2D
  {
    HORZ, VERT
  };

  typedef detran_utils::SP<Sweeper2D>       SP_sweeper;
  typedef State::SP_state                   SP_state;
  typedef detran_utils::InputDB::SP_input   SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef Equation::SP_equation             SP_equation;
  typedef Boundary<_2D>                     Boundary_T;
  typedef Boundary_T::SP_boundary           SP_boundary;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;
  typedef detran_utils::vec_dbl             edge_flux_type;
  typedef Boundary<_2D>::boundary_flux_type boundary_flux_type;

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
            SP_state state,
            SP_boundary boundary,
            SP_sweepsource sweepsource);

  /*!
   *  \brief Sweep over all angles and space.
   *
   */
  virtual inline void sweep(moments_type &phi);

private:

  /// \name Data
  /// \{

  /// Horizontal face fluxes
  edge_flux_type d_psi_h;

  /// Vertical face fluxes
  edge_flux_type d_psi_v;

  int d_side_index[4][2][2];
  /// \}

  /// \name Implementation
  /// \{


  /// \}



};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Sweeper2D.i.hh"

#endif /* SWEEPER2D_HH_ */
