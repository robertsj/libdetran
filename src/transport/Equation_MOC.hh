//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_MOC.hh
 *  @brief  Equation_MOC
 *  @author Jeremy Roberts
 *  @date   Jun 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_MOC_HH_
#define detran_EQUATION_MOC_HH_

#include "transport/transport_export.hh"
#include "DimensionTraits.hh"
#include "angle/ProductQuadrature.hh"
#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "geometry/TrackDB.hh"
#include "geometry/Track.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Equation_MOC
 *  @brief Method of characteristics equation base.
 */
//---------------------------------------------------------------------------//

class TRANSPORT_EXPORT Equation_MOC
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_material::Material::SP_material          SP_material;
  typedef detran_geometry::Mesh::SP_mesh                  SP_mesh;
  typedef detran_geometry::TrackDB::SP_trackdb            SP_trackdb;
  typedef detran_angle::ProductQuadrature::SP_quadrature  SP_quadrature;
  typedef detran_utilities::vec_dbl                       moments_type;
  typedef detran_utilities::vec_dbl                       angular_flux_type;
  typedef detran_utilities::size_t                        size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   */
  Equation_MOC(SP_mesh       mesh,
               SP_material   material,
               SP_quadrature quadrature,
               bool          update_psi)
    :  d_mesh(mesh)
//    ,  d_tracks(mesh->tracks())
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

  virtual ~Equation_MOC(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MOC EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Solve for the cell-center and outgoing edge fluxes.
   *
   *  @param   region      Flat source region (cardinal mesh index)
   *  @param   length      Segment length
   *  @param   width       Track width
   *  @param   source      Reference to sweep source vector for this group
   *  @param   psi_in      Incident flux for this cell
   *  @param   psi_out     Outgoing flux from this cell
   *  @param   phi         Reference to flux moments for this group
   *  @param   psi         Reference to angular flux for this group
   */
  virtual inline void solve(const size_t       region,
                            const double       length,
                            const double       width,
                            moments_type      &source,
                            double            &psi_in,
                            double            &psi_out,
                            moments_type      &phi,
                            angular_flux_type &psi) = 0;

  /**
   *  @brief Setup the equations for a group.
   *  @param g     Current group.
   */
  virtual void setup_group(const size_t g) = 0;

  /**
   *  @brief Setup the equations for an octant.
   *  @param o    Current octant index.
   */
  virtual void setup_octant(const size_t o) = 0;

  /**
   *  @brief Setup the equations for an azimuth.
   *  @param a    Azimuth within octant.
   */
  virtual void setup_azimuth(const size_t a) = 0;

  /**
   *  @brief Setup the equations for a polar angle.
   *  @param p    Polar index.
   */
  virtual void setup_polar(const size_t p) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Problem mesh
  SP_mesh d_mesh;
  /// Tracking data
  SP_trackdb d_tracks;
  /// Material definitions
  SP_material d_material;
  /// MOC Quadrature
  SP_quadrature d_quadrature;
  /// Current mu value
  double d_mu;
  /// current eta value
  double d_eta;
  /// Current ksi value
  double d_xi;
  /// Current track spacing
  double d_spacing;
  /// Inverse of the polar sine
  double d_inv_sin;
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
  /// Current azimuth.
  size_t d_azimuth;
  /// Current polar.
  size_t d_polar;

};


} // end namespace detran

#endif // detran_EQUATION_MOC_HH_

//---------------------------------------------------------------------------//
//              end of file Equation_MOC.hh
//---------------------------------------------------------------------------//
