


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

// Detran
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"
#include "Traits.hh"

// Utilities
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

/*!
 *  \brief Boundary traits for simplify type access.
 */
template <class D>
class EquationTraits
{
public:
  typedef double face_flux_type[D::dimension];
};
template <>
class EquationTraits<_1D>
{
public:
  typedef double face_flux_type;
};


//---------------------------------------------------------------------------//
/*!
 * \class Equation
 * \brief Discrete ordinates equation base.
 */
//---------------------------------------------------------------------------//

template <class D>
class Equation : public Object
{

public:

  typedef SP<Equation>                      SP_equation;
  typedef Material::SP_material             SP_material;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;
  typedef typename
      EquationTraits<D>::face_flux_type     face_flux_type;

  /// Dimension of equation.
  static const int dimension = D::dimension;

  /*!
   *  \brief Constructor
   */
  Equation(SP_mesh mesh,
           SP_material material,
           SP_quadrature quadrature,
           bool update_psi)
    :  d_mesh(mesh)
    ,  d_material(material)
    ,  d_quadrature(quadrature)
    ,  d_update_psi(update_psi)
    ,  d_g(-1)
    ,  d_octant(-1)
    ,  d_angle(-1)
  {
    Require(mesh);
    Require(material);
    Require(quadrature);
    d_mat_map =  mesh->mesh_map("MATERIAL");
    //Ensure(d_mat_map);
  }

  virtual ~Equation(){}

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   \param   i           Cell x index
   *   \param   j           Cell y index
   *   \param   k           Cell z index
   *   \param   source      Cell source
   *   \param   psi_in      Incident flux for this cell
   *   \param   psi_out     Outgoing flux from this cell
   *   \param   phi         Reference to flux moments for this group
   *   \param   psi         Reference to angular flux for this group
   */
  virtual inline void solve(int i,
                            int j,
                            int k,
                            moments_type &source,
                            face_flux_type &psi_in,
                            face_flux_type &psi_out,
                            moments_type &phi,
                            angular_flux_type &psi) = 0;

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
  {return true;}

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
  vec_int d_mat_map;

  /// Update the angular flux?
  bool d_update_psi;

  /// Current group
  int d_g;

  /// Current octant index.
  int d_octant;

  /// Current angle index.
  int d_angle;
};

} // end namespace detran

#endif /* EQUATION_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation.hh
//---------------------------------------------------------------------------//
