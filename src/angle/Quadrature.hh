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

// Other libtran headers
#include "Constants.hh"
#include "Definitions.hh"
#include "DBC.hh"
#include "SP.hh"
#include "Warning.hh"

// System headers
#include <string>

namespace detran
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
class Quadrature : public Object
{

public:

  enum directions
  {
    MU,
    ETA,
    XI,
    END_COSINES
  };

  typedef SP<Quadrature>      SP_quadrature;

  /*!
   *  \brief Constructor.
   *
   *  \param    order       Quadrature order
   *  \param    dim         Spatial dimension
   */
  Quadrature(int order,
             int dim,
             int number_angles,
             std::string name);

  /*!
   *  \brief Pure virtual destructor.
   */
  virtual ~Quadrature() = 0;

  /*!
   *  \brief Return total number of angles.
   */
  int number_angles()
  {
    return d_number_angles;
  }

  /*!
   *  \brief Return total number of octants.
   */
  int number_octants()
  {
    return d_number_octants;
  }

  /*!
   *  \brief Return number of angles per octant.
   */
  int number_angles_octant()
  {
    return d_number_angles_octant;
  }

  /*!
   *  \brief Return cardinal angle index.
   */
  int index(int o, int a)
  {
    Require(o >= 0);
    Require(o < d_number_octants);
    Require(a >= 0);
    Require(a < d_number_angles_octant);
    int angle = a + o * d_number_angles_octant;
    Ensure(angle >= 0);
    Ensure(angle < d_number_angles);
    return angle;
  }

  /*!
   *  \brief Return const reference to weights.
   */
  const vec_dbl& weights()
  {
    return d_weight;
  }

  /*!
   *  \brief Return const reference a cosine vector.
   */
  const vec_dbl& cosines(int dir)
  {
    Require (dir >= 0 and dir < d_dimension);
    if (dir == MU)
      return d_mu;
    else if (dir == ETA)
      return d_eta;
    else
      return d_xi;
  }

  /*!
   *  \brief Return single weight.
   */
  double weight(int a)
  {
    Require(a >= 0);
    Require(a < d_number_angles_octant);
    return d_weight[a];
  }

  /*!
   *  \brief Return single \f$ \mu \f$.
   */
  double mu(int o, int a)
  {
    Require(a >= 0);
    Require(a < d_number_angles_octant);
    Require(o >= 0);
    Require(o < d_number_octants);
    return d_octant_sign[o][MU]*d_mu[a];
  }

  /*!
   *  \brief Return single \f$ \eta \f$.
   */
  double eta(int o, int a)
  {
    Require(a >= 0);
    Require(a < d_number_angles_octant);
    Require(o >= 0);
    Require(o < d_number_octants);
    Require(d_dimension > 1); // 1d calcs have no business with eta...
    return d_octant_sign[o][ETA]*d_eta[a];
  }

  /*!
   *  \brief Return single \f$ \xi \f$.
   */
  double xi(int o, int a)
  {
    Require(a >= 0);
    Require(a < d_number_angles_octant);
    Require(o >= 0);
    Require(o < d_number_octants);
    Require(d_dimension == 3); // 1d/2d calcs have no need for xi...
    return d_octant_sign[o][XI]*d_xi[a];
  }

  /*!
   *  \brief Return one over the integral of unity over all angles.
   */
  static double angular_norm(int d)
  {
    Require(d > 0 and d < 4);
    if (d == 1)
      return 0.5;
    else
      return inv_four_pi;
  }

  /*!
   *  \brief Are the indices valid?
   */
  bool valid_index(int o, int a) const
  {
    if (o >= 0 and o < d_number_octants and
        a >= 0 and a < d_number_angles_octant)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  /*!
   *  \brief Pretty print of the first octant parameters.
   */
  virtual void display() const;

  /*!
   *  \brief Unimplemented DBC method.
   */
  bool is_valid() const
  {
    return true;
  }

protected:

  /// problem dimension
  int d_dimension;

  /// Quadrature order (means different things for each type!)
  int d_order;

  /// number of angles
  int d_number_angles;

  /// number of octants
  int d_number_octants;

  /// number of angle per octant
  int d_number_angles_octant;

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

};

} // end namespace detran

#endif /* QUADRATURE_HH_ */

//---------------------------------------------------------------------------//
//              end of Quadrature.hh
//---------------------------------------------------------------------------//
