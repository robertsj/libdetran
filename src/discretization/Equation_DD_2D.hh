//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_2D.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_2D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_DD_2D_HH_
#define EQUATION_DD_2D_HH_

#include "Equation.hh"

namespace detran
{

/*!
 *  \class Equation_DD_2D
 *  \brief Diamond difference discretization in two dimensions.
 *
 *   See \ref Equation_DD_3D for a general description of the
 *   diamond difference approximation.
 */
class Equation_DD_2D : public Equation<_2D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Equation<_2D> >    SP_equation;
  typedef Equation<_2D>::SP_material              SP_material;
  typedef Equation<_2D>::SP_mesh                  SP_mesh;
  typedef Equation<_2D>::SP_quadrature            SP_quadrature;
  typedef Equation<_2D>::moments_type             moments_type;
  typedef Equation<_2D>::angular_flux_type        angular_flux_type;
  typedef Equation<_2D>::face_flux_type           face_flux_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  Equation_DD_2D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 const bool update_psi);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   See \ref Equation for full description.
   */
  inline void solve(const size_t i,
                    const size_t j,
                    const size_t k,
                    moments_type &source,
                    face_flux_type &psi_in,
                    face_flux_type &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /*!
   *  \brief Setup the equations for a group.
   *  \param g     Current group.
   */
  void setup_group(const size_t g);

  /*!
   *  \brief Setup the equations for an octant.
   *  \param octant    Current octant.
   */
  void setup_octant(const size_t octant);

  /*!
   *  \brief Setup the equations for an angle.
   *  \param angle  Angle index within octant.
   */
  void setup_angle(const size_t angle);

private:

  /// X-directed coefficient, \f$ 2|\mu|/\Delta_x \f$.
  detran_utilities::vec_dbl d_coef_x;

  /// Y-directed coefficient, \f$ 2|\eta|/\Delta_y \f$.
  detran_utilities::vec_dbl d_coef_y;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_DD_2D.i.hh"

#endif /* EQUATION_DD_2D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_2D.hh
//---------------------------------------------------------------------------//
