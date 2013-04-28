//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ProductQuadrature.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  ProductQuadrature class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATURE_HH_
#define detran_angle_PRODUCTQUADRATURE_HH_

#include "Quadrature.hh"

namespace detran_angle
{

/**
 *  @class ProductQuadrature
 *  @brief Base class for product quadratures.
 *
 *  For the discrete ordinates and characteristics methods, a collocation
 *  in angle is used, and approximations of the integral
 *  @f[
 *      \int^{\pi}_0 \sin \phi d\phi \int^{2\pi}_{0} d\theta \,
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
 *  as a product.
 *
 *  As in all quadratures in Detran, we enforce octant symmetry.  Hence,
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
 */
class ANGLE_EXPORT ProductQuadrature: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param dim    Dimension (2 or 3)
   *  @param na     Number of azimuths per octant
   *  @param np     Number of polar angles per octant
   *  @param name   Quadrature name
   */
  ProductQuadrature(const size_t dim,
                    const size_t na,
                    const size_t np,
                    std::string  name);

  /// Pure virtual destructor
  virtual ~ProductQuadrature() = 0;


  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Cardinal angle within octant given azimuth and polar.
  size_t angle(const size_t a, const size_t p) const;
  /// Azimuth index from cardinal within octant
  size_t azimuth(const size_t angle) const;
  /// Polar index from cardinal within octant
  size_t polar(const size_t angle) const;
  double sin_theta(const size_t p) const;
  double cos_theta(const size_t p) const;
  double phi(const size_t a) const;
  double sin_phi(const size_t a) const;
  double cos_phi(const size_t a) const;
  size_t number_azimuths_octant() const;
  size_t number_polar_octant() const;
  double polar_weight(const size_t p) const;
  double azimuth_weight(const size_t a) const;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

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

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Build the product quadrature parameters from polar and azimuth values
  void build_product_quadrature();

private:

  /// Verify the correct ordering.
  void verify() const;

};

} // end namespace detran_angle

#include "ProductQuadrature.i.hh"

#endif /* detran_angle_PRODUCTQUADRATURE_HH_ */
