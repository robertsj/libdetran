//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Quadrature.hh
 * \brief  Quadrature class definition.
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef QUADRATURE_HH_
#define QUADRATURE_HH_

#include "utilities/Constants.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Warning.hh"
#include <string>

/*!
 *  \namespace detran_angle
 *  \brief Contains all angular quadrature and related items for detran
 */
namespace detran_angle
{

//---------------------------------------------------------------------------//
/*!
 *  \class Quadrature
 *  \brief Base quadrature class for discrete ordinates.
 *
 *  All quadratures must be ordered so that the signs of the cosines are
 *  arranged in the following manner:
 *
 *  \verbatim
 *    indices | mu  | eta | xi
 *    -------------------------
 *       1: N | +   |  +  | +    (first octant)
 *     N+1:2N | -   |  +  | +    (second octant)
 *    2N+1:3N | -   |  -  | +    (third octant)
 *    3N+1:4N | +   |  -  | +    (fourth octant)
 *    4N+1:5N | +   |  +  | -    ...
 *    5N+1:6N | -   |  +  | -
 *    6N+1:7N | -   |  -  | -
 *    7N+1:8N | +   |  -  | -
 *  \endverbatim
 *
 *  Note that N is the number of angles per quadrant.  Outside of the given
 *  pattern, the angles need only be consistently ordered, i.e.
 *
 *  \verbatim
 *    abs(mu(i*N+1)) = abs(mu(j*N+1)) for i,j = 0, 1, 2, 3
 *  \endverbatim
 *
 *  though decreasing absolute value is suggested.  This consistency is
 *  helpful for streamlining quadrature sums over all angles.
 *
 *  Because of this ordering, *only the first octant values are stored*.  All
 *  other values are found by multiply by the appropriate octant-dependent
 *  sign.  This assumes, of course, symmetric quadratures.
 *
 */
//---------------------------------------------------------------------------//
class Quadrature
{

public:

  enum directions
  {
    MU,
    ETA,
    XI,
    END_COSINES
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Quadrature>      SP_quadrature;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef detran_utilities::vec2_dbl            vec2_dbl;
  typedef detran_utilities::size_t              size_t;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor.
   *
   *  \param    order           Quadrature order
   *  \param    dim             Spatial dimension
   *  \param    number_angles   Total number of angles
   *  \param    name            Descriptive name
   */
  Quadrature(const size_t order,
             const size_t dim,
             const size_t number_angles,
             const std::string name);

  /// Pure virtual destructor
  virtual ~Quadrature() = 0;

  /// Return total number of angles.
  size_t number_angles() const;

  /// Return total number of octants.
  size_t number_octants() const;

  /// Return number of angles per octant.
  size_t number_angles_octant() const;

  /*!
   *  \brief Return number of azimuths per octant.
   *
   *  This is useful when a product quadrature is defined.
   */
  virtual size_t number_azimuths_octant() const;

  /*!
   *  \brief Return number of polar angles per octant.
   *
   *  This is useful when a product quadrature is defined.
   */
  virtual size_t number_polar_octant() const;

  /*!
   *  \brief Return cardinal angle index.
   *  \param o    Octant index
   *  \param a    Angle within octant
   */
  size_t index(const size_t o, const size_t a);

  /*!
   *  \brief Angle in octant from azimuth and polar
   *
   *  Useful for product quadratures.
   *
   *  \param a  Azimith within octant
   *  \param p  Polar within octant
   */
  virtual size_t angle(const size_t a, const size_t p) const;

  /// Return const reference to weights.
  const vec_dbl& weights() const;

  /*!
   *  \brief Return const reference to a cosine vector
   *  \param dir    Direction of cosine
   */
  const vec_dbl& cosines(const size_t dir) const;

  /*!
   *  \brief Return single weight.
   *  \param a  Angle within octant
   */
  double weight(const size_t a) const;

  /*!
   *  \brief Return single \f$ \mu \f$.
   *  \param o    Octant index
   *  \param a    Angle within octant
   */
  double mu(const size_t o, const size_t a) const;

  /*!
   *  \brief Return single \f$ \eta \f$.
   *  \param o    Octant index
   *  \param a    Angle within octant
   */
  double eta(const size_t o, const size_t a) const;

  /*!
   *  \brief Return single \f$ \xi \f$.
   *  \param o    Octant index
   *  \param a    Angle within octant
   */
  double xi(const size_t o, const size_t a) const;

  /*!
   *  \brief Return one over the integral of unity over all angles.
   *  \param d    Problem dimension
   */
  static double angular_norm(const size_t d);

  /*!
   *  \brief Are the indices valid?
   *  \param o    Octant index
   *  \param a    Angle within octant
   */
  bool valid_index(const size_t o, const size_t a) const;

  size_t dimension() const
  {
    return d_dimension;
  }

  /// Set adjoint.  This changes the octant multipliers.
  void set_adjoint(const bool v = true);

  /// Is adjoint?
  bool is_adjoint() const
  {
    return d_adjoint;
  }

  /*!
   *  \brief Pretty print of the first octant parameters.
   */
  virtual void display() const;

protected:

  /// problem dimension
  size_t d_dimension;

  /// Quadrature order (means different things for each type!)
  size_t d_order;

  /// number of angles
  size_t d_number_angles;

  /// number of octants
  size_t d_number_octants;

  /// number of angle per octant
  size_t d_number_angles_octant;

  /// Quadrature weight
  vec_dbl d_weight;

  /// x-axis cosines
  vec_dbl d_mu;

  /// y-axis cosines
  vec_dbl d_eta;

  /// z-axis cosines
  vec_dbl d_xi;

  /// type of quadrature being used
  std::string d_name;

  /// octant cosign signs (e.g [1 1 1] for octant 1)
  vec2_dbl d_octant_sign;

  /// Is this an adjoint problem?
  bool d_adjoint;

};

} // end namespace detran_angle

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//
#include "Quadrature.i.hh"

#endif /* QUADRATURE_HH_ */

//---------------------------------------------------------------------------//
//              end of Quadrature.hh
//---------------------------------------------------------------------------//
