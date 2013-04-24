//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_HH_
#define detran_EQUATION_HH_

#include "transport/transport_export.hh"
#include "DimensionTraits.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "angle/Quadrature.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

/**
 *  @namespace  detran
 *  @brief      Core elements for discretizing and solving transport equations
 */
namespace detran
{

/**
 *  @class Equation traits
 *  @brief Traits for defining the face flux type for a discretization
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
/**
 *  @class Equation
 *  @brief Discrete ordinates equation base.
 *  @tparam  D   Problem dimension
 */
//---------------------------------------------------------------------------//
template <class D>
class Equation
{

public:

  /// Dimension of equation.
  static const int dimension = D::dimension;

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Equation>                  SP_equation;
  typedef detran_material::Material::SP_material          SP_material;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature         SP_quadrature;
  typedef detran_utilities::vec_dbl                       moments_type;
  typedef detran_utilities::vec_dbl                       angular_flux_type;
  typedef typename EquationTraits<D>::face_flux_type      face_flux_type;
  typedef detran_utilities::size_t                        size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mesh           Geometry
   *  @param material       Material database
   *  @param quadrature     Angular mesh
   *  @param update_psi     Flag for keeping the angular flux
   */
  Equation(SP_mesh mesh,
           SP_material material,
           SP_quadrature quadrature,
           const bool update_psi)
    :  d_mesh(mesh)
    ,  d_material(material)
    ,  d_quadrature(quadrature)
    ,  d_update_psi(update_psi)
    ,  d_g(0)
    ,  d_octant(0)
    ,  d_angle(0)
  {
    // Preconditions
    Require(mesh);
    Require(material);
    Require(quadrature);

    // Get the material map.
    d_mat_map =  mesh->mesh_map("MATERIAL");
  }

  // Virtual destructor
  virtual ~Equation(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *   @brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   @param   i           Cell x index
   *   @param   j           Cell y index
   *   @param   k           Cell z index
   *   @param   source      Reference to source vector
   *   @param   psi_in      Incident flux for this cell
   *   @param   psi_out     Outgoing flux from this cell
   *   @param   phi         Reference to flux moments for this group
   *   @param   psi         Reference to angular flux for this group
   */
  virtual void solve(const size_t i,
                     const size_t j,
                     const size_t k,
                     moments_type &source,
                     face_flux_type &psi_in,
                     face_flux_type &psi_out,
                     moments_type &phi,
                     angular_flux_type &psi) = 0;

  /**
   *  @brief Setup the equations for a group.
   *  @param g     Current group.
   */
  virtual void setup_group(const size_t g) = 0;

  /**
   *  @brief Setup the equations for an octant.
   *  @param octant    Current octant.
   */
  virtual void setup_octant(const size_t octant) = 0;

  /**
   *  @brief Setup the equations for an angle.
   *  @param angle  Angle index within octant
   */
  virtual void setup_angle(const size_t angle) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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
  detran_utilities::vec_int d_mat_map;
  /// Update the angular flux?
  bool d_update_psi;
  /// Current group
  size_t d_g;
  /// Current octant index.
  size_t d_octant;
  /// Current angle index.
  size_t d_angle;

};

} // end namespace detran

#endif /* detran_EQUATION_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation.hh
//---------------------------------------------------------------------------//
