/*
 * PolarQuadrature.hh
 *
 *  Created on: Jun 22, 2012
 *      Author: robertsj
 */

#ifndef POLARQUADRATURE_HH_
#define POLARQUADRATURE_HH_

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

// System
#include <vector>

namespace detran
{

class PolarQuadrature: public Object
{

public:

  typedef SP<PolarQuadrature>   SP_polar;

  /*!
   *  \brief Constructor
   *  \param number_polar   Number of angles in positive half space
   */
  PolarQuadrature(int number_polar)
    : d_number_polar(number_polar)
    , d_sin(number_polar, 0.0)
    , d_cos(number_polar, 0.0)
    , d_weight(number_polar, 0.0)
  {
    Require(number_polar > 0);
  }

  /// Virtual destructor.
  virtual ~PolarQuadrature(){};

  int number_polar() const
  {
    return d_number_polar;
  }

  /*!
   *  \brief Return a polar sine
   *  \param a
   */
  double sin_theta(int a) const
  {
    Require(a < d_number_polar);
    return d_sin[a];
  }

  /*!
   *  \brief Return a polar cosine
   *  \param a
   */
  double cos_theta(int a) const
  {
    Require(a < d_number_polar);
    return d_cos[a];
  }

  /*!
   *  \brief Return a polar weight
   *  \param a
   */
  double weight(int a) const
  {
    Require(a < d_number_polar);
    return d_weight[a];
  }

  /// DBC function
  bool is_valid() const
  {
    return true;
  }

protected:

  /// \name Protected Data
  /// \{

  /// Number of polar angles for a single half space.
  int d_number_polar;

  /// Vector of sines
  std::vector<double> d_sin;

  /// Vector of cosines
  std::vector<double> d_cos;

  /// Vector of weights (sums to 1).
  std::vector<double> d_weight;

  /// \}

};

}


#endif /* POLARQUADRATURE_HH_ */
