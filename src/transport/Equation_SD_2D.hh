//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SD_2D.hh
 * \author robertsj
 * \date   Jun 8, 2012
 * \brief  Equation_SD_2D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SD_2D_HH_
#define EQUATION_SD_2D_HH_

// Detran headers
#include "Equation.hh"

// Detran utilities
#include "Definitions.hh"

namespace detran
{

/*!
 *  \class Equation_SD_2D
 *  \brief Step difference discretization in two dimensions.
 *
 */
class Equation_SD_2D : public Equation<_2D>
{

public:

  typedef SP<Equation<_2D> >                SP_equation;
  typedef Equation<_2D>::SP_material        SP_material;
  typedef Equation<_2D>::SP_mesh            SP_mesh;
  typedef Equation<_2D>::SP_quadrature      SP_quadrature;
  typedef Equation<_2D>::moments_type       moments_type;
  typedef Equation<_2D>::angular_flux_type  angular_flux_type;
  typedef Equation<_2D>::face_flux_type     face_flux_type;

  /*!
   *  \brief Constructor
   */
  Equation_SD_2D(SP_mesh mesh,
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

  /// Y-directed coefficient, \f$ 2|\eta|/\Delta_y \f$.
  vec_dbl d_coef_y;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SD_2D.i.hh"

#endif /* EQUATION_SD_2D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SD_2D.hh
//---------------------------------------------------------------------------//



