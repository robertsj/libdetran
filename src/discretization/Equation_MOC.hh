//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_MOC.hh
 * \brief  Equation_MOC 
 * \author Jeremy Roberts
 * \date   Jun 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_MOC_HH_
#define EQUATION_MOC_HH_

// Detran
#include "Material.hh"
#include "MeshMOC.hh"
#include "QuadratureMOC.hh"
#include "State.hh"
#include "TrackDB.hh"
#include "Track.hh"
#include "Traits.hh"

// Utilities
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Equation_MOC
 * \brief Method of characteristics equation base.
 */
//---------------------------------------------------------------------------//

class Equation_MOC
{

public:

  typedef Material::SP_material             SP_material;
  typedef MeshMOC::SP_mesh                  SP_mesh;
  typedef TrackDB::SP_trackdb               SP_trackdb;
  typedef QuadratureMOC::SP_quadrature      SP_quadrature;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

  /*!
   *  \brief Constructor
   */
  Equation_MOC(SP_mesh mesh,
               SP_material material,
               SP_quadrature quadrature,
               bool update_psi)
    :  d_mesh(mesh)
    ,  d_tracks(mesh->tracks())
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

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   \param   region      Flat source region (cardinal mesh index)
   *   \param   length      Segment length
   *   \param   source      Reference to sweep source vector for this group
   *   \param   psi_in      Incident flux for this cell
   *   \param   psi_out     Outgoing flux from this cell
   *   \param   phi         Reference to flux moments for this group
   *   \param   psi         Reference to angular flux for this group
   */
  virtual inline void solve(int region,
                            double length,
                            moments_type &source,
                            double &psi_in,
                            double &psi_out,
                            moments_type &phi,
                            angular_flux_type &psi) = 0;


  /*!
   *  \brief Setup the equations for a group.
   *  \param g     Current group.
   */
  virtual void setup_group(int g) = 0;

  /*!
   *  \brief Setup the equations for an octant.
   *  \param o    Current octant index.
   */
  virtual void setup_octant(int o) = 0;

  /*!
   *  \brief Setup the equations for an azimuth.
   *  \param a    Azimuth within octant.
   */
  virtual void setup_azimuth(int a) = 0;

  /*!
   *  \brief Setup the equations for a polar angle.
   *  \param p    Polar index.
   */
  virtual void setup_polar(int p) = 0;

  /// \}

  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  /// \{

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
  vec_int d_mat_map;

  /// Update the angular flux?
  bool d_update_psi;

  /// Current group
  int d_g;

  /// Current octant index.
  int d_octant;

  /// Current angle index.
  int d_angle;

  /// Current azimuth.
  int d_azimuth;

  /// Current polar.
  int d_polar;

  /// \}

};


} // end namespace detran

#endif // EQUATION_MOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file Equation_MOC.hh
//---------------------------------------------------------------------------//
