//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_HH_
#define EQUATION_HH_

// Detran utilities
#include "Definitions.hh"
#include "SP.hh"
//#include "InputDB.hh"

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Equation
 * \brief Discrete ordinates equation base.
 */
//---------------------------------------------------------------------------//
class Equation : public detran_utils::Object
{

public:

  typedef detran_utils::SP<Equation>        SP_equation;
  typedef Material::SP_material             SP_material;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;
  typedef double                            face_flux_1d[1];
  typedef double                            face_flux_2d[2];
  typedef double                            face_flux_3d[3];


  /*!
   *  \brief Constructor
   */
  Equation(SP_mesh mesh,
           SP_material material,
           SP_quadrature quadrature,
           bool store_psi = false)
    :  d_mesh(mesh)
    ,  d_material(material)
    ,  d_quadrature(quadrature)
    ,  d_store_psi(store_psi)
  {
    Require(mesh);
    Require(material);
    Require(quadrature);
    d_mat_map =  mesh->mesh_map("MATERIAL");
    //Ensure(d_mat_map);
  }

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   Note, updating psi is optional.  Only if \ref store_psi is
   *   true will it be accesses; if store_psi is false, it's
   *   assumed the corresponding vectors are not sized.
   *
   *   This would be virtual, but template functions cannot be.  Hence,
   *   we simply make it protected so that public access is via
   *   derived implementations.
   *
   *   \param   g           Group
   *   \param   i           Cell x index
   *   \param   j           Cell y index
   *   \param   k           Cell z index
   *   \param   psi_in      Incident flux for this cell
   *   \param   source      Cell source
   *   \param   psi_out     Outgoing flux from this cell
   *   \param   phi         Reference to flux moments for this group
   *   \param   psi         Reference to angular flux for this group
   *   \tparam  F           Edge flux type.  A vector of 1, 2, or 3
   *                        doubles depending on dimension.
   */
  template <class F>
  void solve(int g,
             int i,
             int j,
             int k,
             F psi_in,
             double source,
             F psi_out,
             moments_type &phi,
             angular_flux_type &psi)
  {
    /* ... */
  }

public:

  /*!
   *  @brief Setup the equations for a group.
   *  @param g     Current group.
   */
  virtual void setup_group(int g) = 0;

  /*!
   *  @brief Setup the equations for an octant.
   *  @param octant    Current octant.
   */
  virtual void setup_octant(int octant) = 0;

  /*!
   *  \brief Setup the equations for an angle.
   *  \param angle  Angle index within octant
   */
  virtual void setup_angle(int angle) = 0;

  /// \}

  bool is_valid() const
  { /* ... */ }

protected:

  /// Problem mesh
  SP_mesh d_mesh;

  /// Material definitions
  SP_material d_material;

  /// Quadrature
  SP_quadrature d_quadrature;

  /// Current mu value
  double d_mu;

  /// current eta value
  double d_eta;

  /// Current ksi value
  double d_ksi;

  /// Material map
  detran_utils::vec_int d_mat_map;

  /// Update the angular flux?
  bool d_store_psi;
};

} // end namespace detran

#endif /* EQUATION_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation.hh
//---------------------------------------------------------------------------//
