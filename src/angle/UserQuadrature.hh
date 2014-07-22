//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  UserQuadrature.hh
 *  @brief UserQuadrature class definition
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_UserQuadrature_HH_
#define detran_angle_UserQuadrature_HH_

#include "Quadrature.hh"
#include <cmath>

namespace detran_angle
{

/**
 *  @class UserQuadrature
 *  @brief 1D/2D/3D quadrature based on user-defined abscissa and weights
 */
class ANGLE_EXPORT UserQuadrature : public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param np     Number of abscissa per axis in one octant (2*np = order)
   *  @param dim    Problem dimension
   */
  UserQuadrature(size_t         dim,
                 const vec_dbl &mu_v,
                 const vec_dbl &eta_v,
                 const vec_dbl &wt_v)
    : Quadrature(dim, mu_v.size() * std::pow(2, dim), "user")
  {
    Require( (dim == 1) or ((dim > 1) & (mu_v.size() == eta_v.size())) );
    Require(mu_v.size() > 0);
    Require(mu_v.size() == wt_v.size());
    for (int i = 0; i < d_number_angles_octant; ++i)
    {
      std::cout << mu_v[i] << " " << eta_v[i] << std::endl;
      d_mu[i] = mu_v[i];
      if (dim > 1)
      {
        d_eta[i] = eta_v[i];
        d_xi[i]  = std::sqrt(1.0 - d_mu[i]*d_mu[i] - d_eta[i]*d_eta[i]);
      }
      d_weight[i] = wt_v[i];
    }
  }

private:

};



} // end namespace detran_angle

#endif /* detran_angle_UserQuadrature_HH_ */

//---------------------------------------------------------------------------//
//              end of Level_Symmetric_3D.hh
//---------------------------------------------------------------------------//
