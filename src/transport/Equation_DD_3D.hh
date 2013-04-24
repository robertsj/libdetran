//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_3D.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_3D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_DD_3D_HH_
#define detran_EQUATION_DD_3D_HH_

#include "Equation.hh"

namespace detran
{

/**
 *  @class Equation_DD_3D
 *  @brief Diamond difference discretization in three dimensions.
 *
 *  Fill me.
 *
 */
class TRANSPORT_EXPORT Equation_DD_3D : public Equation<_3D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Equation<_3D> >    SP_equation;
  typedef Equation<_3D>::SP_material              SP_material;
  typedef Equation<_3D>::SP_mesh                  SP_mesh;
  typedef Equation<_3D>::SP_quadrature            SP_quadrature;
  typedef Equation<_3D>::moments_type             moments_type;
  typedef Equation<_3D>::angular_flux_type        angular_flux_type;
  typedef Equation<_3D>::face_flux_type           face_flux_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  Equation_DD_3D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 bool update_psi);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Solve for the cell-center and outgoing edge fluxes.
  inline void solve(const size_t i,
                    const size_t j,
                    const size_t k,
                    moments_type &source,
                    face_flux_type &psi_in,
                    face_flux_type &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /// Setup the equations for a group.
  void setup_group(const size_t g);

  /// Setup the equations for an octant.
  void setup_octant(const size_t octant);

  /// Setup the equations for an angle.
  void setup_angle(const size_t angle);


private:

  /// X-directed coefficient, \f$ 2|\mu|/\Delta_x \f$.
  detran_utilities::vec_dbl d_coef_x;

  /// Y-directed coefficient, \f$ 2|\eta|/\Delta_y \f$.
  detran_utilities::vec_dbl d_coef_y;

  /// Z-directed coefficient, \f$ 2|\xi|/\Delta_z \f$.
  detran_utilities::vec_dbl d_coef_z;
};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_DD_3D.i.hh"

#endif /* detran_EQUATION_DD_3D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_3D.hh
//---------------------------------------------------------------------------//
