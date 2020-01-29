//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Quadrature.hh
 *  @brief Quadrature class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_QUADRATURE_HH_
#define detran_angle_QUADRATURE_HH_

#include "angle/angle_export.hh"
#include "utilities/Constants.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include "utilities/Warning.hh"
#include <string>

/**
 *  @namespace detran_angle
 *  @brief Contains all angular quadrature and related items for detran
 */
namespace detran_angle
{

/**
 *  @class Quadrature
 *  @brief Base quadrature class for transport calculations
 *
 *  For all Detran quadratures, octant symmetry is assumed.  As a result,
 *  only the first octant abscissa and weights are stored.
 *
 *  The octants are ordered as follows, with
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
 *  other values are found by multiplying by the appropriate octant-dependent
 *  sign.
 *
 */

class ANGLE_EXPORT Quadrature
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  enum directions
  {
    MU, ETA, XI, END_COSINES
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Quadrature>      SP_quadrature;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef detran_utilities::vec2_dbl            vec2_dbl;
  typedef detran_utilities::vec_int             vec_int;
  typedef detran_utilities::vec2_int            vec2_int;
  typedef detran_utilities::size_t              size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param    dim             Spatial dimension
   *   @param    number_angles   Total number of angles
   *   @param    name            Descriptive name
   */
  Quadrature(const size_t      dim,
             const size_t      number_angles,
             const std::string name);

  /// Pure virtual destructor
  virtual ~Quadrature() = 0;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Return total number of angles.
  size_t number_angles() const;
  /// Return total number of octants.
  size_t number_octants() const;
  /// Return number of angles per octant.
  size_t number_angles_octant() const;

  /**
   *  @brief Return cardinal angle index.
   *  @param o    Octant index
   *  @param a    Angle within octant
   */
  size_t index(const size_t o, const size_t a);

  /// Return const reference to weights.
  const vec_dbl& weights() const;

  /**
   *  @brief Return const reference to a cosine vector
   *  @param dir    Direction of cosine
   */
  const vec_dbl& cosines(const size_t dir) const;

  /**
   *  @brief Return single weight.
   *  @param a  Angle within octant
   */
  double weight(const size_t a) const;

  /**
   *  @brief Return single \f$ \mu \f$.
   *  @param o    Octant index
   *  @param a    Angle within octant
   */
  double mu(const size_t o, const size_t a) const;

  /**
   *  @brief Return single \f$ \eta \f$.
   *  @param o    Octant index
   *  @param a    Angle within octant
   */
  double eta(const size_t o, const size_t a) const;

  /**
   *  @brief Return single \f$ \xi \f$.
   *  @param o    Octant index
   *  @param a    Angle within octant
   */
  double xi(const size_t o, const size_t a) const;

  /**
   *  @brief Return one over the integral of unity over all angles.
   *  @param d    Problem dimension
   */
  static double angular_norm(const size_t d);

  /// Get vector of incident octants for a side
  const vec_int& incident_octant(const size_t s);

  /// Get vector of outgoing octants for a side
  const vec_int& outgoing_octant(const size_t s);

  /**
   *  @brief Are the indices valid?
   *  @param o    Octant index
   *  @param a    Angle within octant
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

  /// Pretty print of the first octant parameters.
  virtual void display() const;


protected:

  /// problem dimension
  size_t d_dimension;
  /// number of angles
  size_t d_number_angles;
  /// number of octants
  size_t d_number_octants;
  /// number of angle per octant
  size_t d_number_angles_octant;
  /// quadrature weight
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
  /// incident octants
  vec2_int d_incident_octants;
  /// outgoing octants
  vec2_int d_outgoing_octants;
  /// Is this an adjoint problem?
  bool d_adjoint;

};

ANGLE_TEMPLATE_EXPORT(detran_utilities::SP<Quadrature>)

} // end namespace detran_angle

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "Quadrature.i.hh"

#endif /* detran_angle_QUADRATURE_HH_ */

//----------------------------------------------------------------------------//
//              end of Quadrature.hh
//----------------------------------------------------------------------------//
