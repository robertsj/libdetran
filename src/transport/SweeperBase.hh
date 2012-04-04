//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SweeperBase.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  SweeperBase class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPERBASE_HH_
#define SWEEPERBASE_HH_

// Other libtran headers
#include "Equation.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"

// Other libtran utils headers
#include "Definitions.hh"
#include "SP.hh"
#include "InputDB.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class SweeperBase
 * \brief Base sweeper for discrete ordinates problems.
 *
 * The within-group transport equation is
 * \f[
 *      \mathbf{L}\psi = Q \, ,
 * \f]
 * where \f$ \mathbf{L} \f$ is the streaming and collision operator and
 * \f$ Q \f$ is a discrete representation of all source contributions.
 *
 * To invert the operator \f$ \mathbf{L} \f$, we "sweep" over the mesh for all
 * angles,
 * which gives us updated angular fluxes in each cell.  Actually, the
 * flux *moments* are updated, while the discrete angular flux is
 * optionally stored.
 *
 */
//---------------------------------------------------------------------------//
class SweeperBase : public detran_utils::Object
{

public:

  typedef detran_utils::SP<SweeperBase>     SP_sweeper;
  typedef State::SP_state                   SP_state;
  typedef detran_utils::InputDB::SP_input   SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef Equation::SP_equation             SP_equation;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   */
  SweeperBase(SP_input input,
              SP_mesh mesh,
              SP_material material,
              SP_quadrature quadrature,
              SP_state state)
    :  d_input(input)
    ,  d_mesh(mesh)
    ,  d_quadrature(quadrature)
    ,  d_state(state)
    ,  d_g(-1)
    ,  d_update_psi(false)
  {
    Require(input);
    Require(mesh);
    Require(quadrature);
    Require(state);

    std::string equation;
    if (input->check("update_psi"))
    {
      d_update_psi = input->get<int>("update_psi")==1;
    }
  }

  /*!
   *  \brief Sweep over all angles and space.
   *
   *  Note, if the angular flux is to be updated,
   *  it is done directly to via State.  Having
   *  sweep take the flux as an explicit argument
   *  allows various input types (e.g. Krylov
   *  vectors) without having to go through State.
   *
   *  Note, here the source is limited to a moments-based
   *  one.  We may want a discrete source later; that can
   *  be added with a third argument of psi type.
   *
   */
  virtual inline void sweep(moments_type &phi,
                            moments_type &source) = 0;

  /// Setup the equations for the group
  void setup_group(int g)
  {
    d_g = g;
    d_equation->setup_group(g);
  }

  /// Allows the psi update to occur whenever needed
  void set_update_psi(bool v)
  {
    d_update_psi = v;
  }

  virtual bool is_valid() const
  {
    /* ... */
  }

protected:

  /// Input database
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// State vectors
  SP_state d_state;

  /// Equation
  SP_equation d_equation;

  /// Current group
  int d_g;

  /// Update the angular flux?
  bool d_update_psi;

};

} // end namespace detran

#endif /* SWEEPERBASE_HH_ */
