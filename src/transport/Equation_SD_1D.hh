//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SD_1D.hh
 * \author robertsj
 * \date   Jun 9, 2012
 * \brief  Equation_SD_1D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SD_1D_HH_
#define EQUATION_SD_1D_HH_

// Detran headers
#include "Equation.hh"

// Detran utilities
#include "Definitions.hh"

namespace detran
{

/*!
 *  \class Equation_SD_1D
 *  \brief Step difference discretization in one dimension.
 *
 *  The step difference approximation defines
 *  \f[
 *      \psi_{i,n} = \left\{
 *        \begin{array}{l l}
 *          \psi_{i+1/2,n}      & \quad \text{if $\mu_n > 0$} \\
 *          \psi_{i-1/2,n}      & \quad \text{if $\mu_n < 0$} \\
 *        \end{array} \right\}
 *  \f]
 *
 *  This is a first order method, but is strictly positive.
 *  Moreover, it essentially defines a "flat flux" within
 *  the cell, and so it is consistent for DGM.
 *
 */
class Equation_SD_1D : public Equation<_1D>
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
  Equation_SD_1D(SP_mesh mesh,
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
   *  \brief Setup the equations for a group.
   *  \param g     Current group.
   */
  void setup_group(int g);

  /*!
   *  \brief Setup the equations for an octant.
   *  \param octant    Current octant.
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

#include "Equation_SD_1D.i.hh"

#endif /* EQUATION_SD_1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SD_1D.hh
//---------------------------------------------------------------------------//
