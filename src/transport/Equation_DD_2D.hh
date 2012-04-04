//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_DD_2D.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_DD_2D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_DD_2D_HH_
#define EQUATION_DD_2D_HH_

// Detran headers
#include "Equation.hh"

// Detran utilities
#include "Definitions.hh"

namespace detran
{

class Equation_DD_2D : public Equation
{

public:

  typedef detran_utils::SP<Equation>        SP_equation;
  typedef Equation::SP_material             SP_material;
  typedef Equation::SP_mesh                 SP_mesh;
  typedef Equation::SP_quadrature           SP_quadrature;
  typedef Equation::moments_type            moments_type;
  typedef Equation::angular_flux_type       angular_flux_type;
  typedef Equation::face_flux_2d            face_flux_type;

  /*!
   *  \brief Constructor
   */
  Equation_DD_2D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 bool update_psi);

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   Note, updating psi is optional.  Only if \ref store_psi is
   *   true will it be accesses; if store_psi is false, it's
   *   assumed the corresponding vectors are not sized.
   *
   *   Also, the
   *
   *   \param   g           Group
   *   \param   i           Cell x index
   *   \param   j           Cell y index
   *   \param   k           Cell z index
   *   \param   source      Cell source
   *   \param   psi_in      Incident flux for this cell
   *   \param   psi_out     Outgoing flux from this cell
   *   \param   phi         Reference to flux moments for this group
   *   \param   psi         Reference to angular flux for this group
   *   \tparam  F           Edge flux type.  A vector of 1, 2, or 3
   *                        doubles depending on dimension.
   */
  template <class F>
  inline void solve(int i,
                    int j,
                    int k,
                    moments_type &source,
                    F &psi_in,
                    F &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /*!
   *  @brief Setup the equations for a group.
   *  @param g     Current group.
   */
  virtual void setup_group(int g);

  /*!
   *  @brief Setup the equations for an octant.
   *  @param octant    Current octant.
   */
  virtual void setup_octant(int octant);

  /*!
   *  \brief Setup the equations for an angle.
   *  \param angle  Angle index within octant.
   */
  virtual void setup_angle(int angle);

  /// \}

private:

  /// X-directed coefficient, \f$ 2|\mu|/\Delta_x \f$.
  detran_utils::vec_dbl d_coef_x;

  /// Y-directed coefficient, \f$ 2|\eta|/\Delta_y \f$.
  detran_utils::vec_dbl d_coef_y;

  /// Current octant index.
  int d_octant;

  /// Current angle index.
  int d_angle;

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
