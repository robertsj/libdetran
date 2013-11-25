//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.hh
 *  @brief ProductQuadrature class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATURE_HH_
#define detran_angle_PRODUCTQUADRATURE_HH_

#include "angle/Quadrature.hh"
#include "angle/PolarQuadrature.hh"
#include "angle/AzimuthalQuadrature.hh"

namespace detran_angle
{

/**
 *  @class ProductQuadrature
 *  @brief Base class for product quadratures
 *
 *  For the discrete ordinates and characteristics methods, a collocation
 *  in angle is used, and approximations of the integral
 *  @f[
 *      \int^{\pi}_0 \sin \theta d\theta \int^{2\pi}_{0} d\phi \,
 *  @f]
 *  are needed over a discrete set of angular abscissa @f$ \hat{\Omega}_n @f$
 *  where
 *  @f[
 *      \hat{\Omega}_n = \hat{i} \sin \phi_n \cos \theta_n +
 *                       \hat{j} \sin \phi_n \sin \theta_n +
 *                       \hat{k} \cos \phi_n
 *  @f]
 *  with associated weights @f$ w_n @f$.
 *
 *  A product quadrature separates the angular variables, defining
 *  one dimensional quadrature formulas for each and combining the two
 *  as an outer product.
 *
 *  As in all quadratures in Detran, octant symmetry is enforced.  Hence,
 *  the azimuthal and polar quadratures are defined only within an
 *  octant.
 *
 *  Furthermore, product quadratures must order angles so that
 *  the polar angle changes most quickly and in increasing order.
 *  For a given polar level, the azimuth (from x) increases from
 *  0 to pi/2,  i.e. mu goes from 1 to 0.  This is enforced so that
 *  applications dependent on ordering need not be responsible.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant    -- number polar angles per octant
 *    - quad_number_azimuth_octant  -- number azimuths per octant
 *
 *  @tparam   A     Azimuthal Quadrature
 *  @tparam   P     Polar Quadrature
 */
class ANGLE_EXPORT ProductQuadrature: public Quadrature
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<ProductQuadrature>   SP_quadrature;

  struct index_t
  {
    size_t octant;
    size_t angle;
    size_t io_octant;
  };

  typedef std::vector<index_t>          vec1_index_t;
  typedef std::vector<vec1_index_t>     vec2_index_t;
  typedef std::vector<vec2_index_t>     vec3_index_t;
  typedef detran_utilities::vec_size_t  vec_size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    dim       Dimension (2 or 3)
   *  @param    na        Number of azimuths per octant
   *  @param    np        Number of polar angles per octant
   *  @param    name      Quadrature name
   *  @param    normalize Ensure normalization of polar and azimuthal sets
   */
  ProductQuadrature(const size_t       dim,
                    const size_t       na,
                    const size_t       np,
                    const std::string &name = "product",
                    const bool         normalize = false);

  /// Virtual destructor
  virtual ~ProductQuadrature() = 0;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  //@{
  ///  Azimuthal angle, cosine, sine, or weight via azimuthal index
  double phi(const size_t a) const;
  double cos_phi(const size_t a) const;
  double sin_phi(const size_t a) const;
  double azimuth_weight(const size_t a) const;
  //@}

  //@{
  ///  Polar cosine, sine  or weight via polar index
  double cos_theta(const size_t p) const;
  double sin_theta(const size_t p) const;
  double polar_weight(const size_t p) const;
  //@}

  /// Number of azimuths per octant
  size_t number_azimuths_octant() const {return d_number_azimuth_octant;}
  /// Number of polar angles per octant
  size_t number_polar_octant() const {return d_number_polar_octant;}
  /// Number of incident or outgoing azimuths along a side
  size_t number_azimuths(const size_t s) const;
  /// Number of incident or outgoing polar angles along a side
  size_t number_polar(const size_t s) const;

  /// Cardinal angle within octant given azimuth and polar.
  size_t angle(const size_t a, const size_t p) const;
  /// Azimuth index from cardinal within octant
  size_t azimuth(const size_t angle) const;
  /// Polar index from cardinal within octant
  size_t polar(const size_t angle) const;

  // Prepend class on the following because onld SWIG fails on this
  // particular nested type.

  /// Incident octant/angle given cardinal azimuth and polar index on a side
  ProductQuadrature::index_t
  incident_index(const size_t s, const size_t a, const size_t p) const;
  /// Outgoing octant/angle given cardinal azimuth and polar index on a side
  ProductQuadrature::index_t
  outgoing_index(const size_t s, const size_t a, const size_t p) const;


  void display_indices();

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Number of azimuths per octant
  size_t d_number_azimuth_octant;
  /// Number of polar angles per octant
  size_t d_number_polar_octant;
  /// Azimuths
  vec_dbl d_phi;
  /// Cosine of azimuths
  vec_dbl d_cos_phi;
  /// Sine of azimuths
  vec_dbl d_sin_phi;
  /// Azimuth weights
  vec_dbl d_azimuth_weight;
  /// Cosine of polar angles
  vec_dbl d_cos_theta;
  /// Sine of polar angles
  vec_dbl d_sin_theta;
  /// Polar weights
  vec_dbl d_polar_weight;
  ///
  vec_size_t d_number_polar;
  vec_size_t d_number_azimuths;
  /// Angle indices for each side [side][azimuth][polar]
  vec3_index_t d_incident_indices;
  vec3_index_t d_outgoing_indices;


  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Build azimuth and polar sets with optional normalization
  void build();
  /// Verify the correct ordering.
  void verify() const;

};

} // end namespace detran_angle

#endif /* detran_angle_PRODUCTQUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.hh
//----------------------------------------------------------------------------//
