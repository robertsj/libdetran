//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_1D.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_1D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_DD_1D_HH_
#define EQUATION_DD_1D_HH_

// Detran headers
#include "Equation.hh"

// Detran utilities
#include "Definitions.hh"

namespace detran
{

/*!
 *  \class Equation_DD_1D
 *  \brief Diamond difference discretization in two dimensions.
 *
 */
class Equation_DD_1D : public Equation<_1D>
{

public:

  typedef SP<Equation<_1D> >                SP_equation;
  typedef Equation<_1D>::SP_material        SP_material;
  typedef Equation<_1D>::SP_mesh            SP_mesh;
  typedef Equation<_1D>::SP_quadrature      SP_quadrature;
  typedef Equation<_1D>::moments_type       moments_type;
  typedef Equation<_1D>::angular_flux_type  angular_flux_type;
  typedef Equation<_1D>::face_flux_type     face_flux_type;

  /*!
   *  \brief Constructor
   */
  Equation_DD_1D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 bool update_psi);

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   See \ref Equation for full description.
   */
  inline void solve(int i,
                    int j,
                    int k,
                    moments_type &source,
                    face_flux_type &psi_in,
                    face_flux_type &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /*!
   *  @brief Setup the equations for a group.
   *  @param g     Current group.
   */
  void setup_group(int g);

  /*!
   *  @brief Setup the equations for an octant.
   *  @param octant    Current octant.
   */
  void setup_octant(int octant);

  /*!
   *  \brief Setup the equations for an angle.
   *  \param angle  Angle index within octant.
   */
  void setup_angle(int angle);

  /// \}

private:

  /// X-directed coefficient, \f$ 2|\mu|/\Delta_x \f$.
  vec_dbl d_coef_x;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_DD_1D.i.hh"

#endif /* EQUATION_DD_1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_1D.hh
//---------------------------------------------------------------------------//
