//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  ProductQuadratureAP.hh
 *  @brief ProductQuadratureAP class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATUREAP_HH_
#define detran_angle_PRODUCTQUADRATUREAP_HH_

#include "ProductQuadrature.hh"

namespace detran_angle
{

/**
 *  @class ProductQuadratureAP
 *  @brief Product quadrature templated on azimuthal and polar quadratures
 */

template <class A, class P>
class ANGLE_EXPORT ProductQuadratureAP: public ProductQuadrature
{

public:

  /// Constructor
  ProductQuadratureAP(const size_t dim,
                      const size_t na,
                      const size_t np,
                      const bool   normalize = false);

  ~ProductQuadratureAP(){/* ... */}

};

//----------------------------------------------------------------------------//
// CONVENIENCE TYPEDEFS
//----------------------------------------------------------------------------//

typedef ProductQuadratureAP<AzimuthalGL, PolarGL>         Product_GL_GL;
typedef ProductQuadratureAP<AzimuthalGL, PolarDGL>        Product_GL_DGL;
typedef ProductQuadratureAP<AzimuthalGL, PolarU>          Product_GL_U;
typedef ProductQuadratureAP<AzimuthalGL, PolarTY>         Product_GL_TY;
//
typedef ProductQuadratureAP<AzimuthalDGL, PolarGL>        Product_DGL_GL;
typedef ProductQuadratureAP<AzimuthalDGL, PolarDGL>       Product_DGL_DGL;
typedef ProductQuadratureAP<AzimuthalDGL, PolarU>         Product_DGL_U;
typedef ProductQuadratureAP<AzimuthalDGL, PolarTY>        Product_DGL_TY;
//
typedef ProductQuadratureAP<AzimuthalU, PolarGL>          Product_U_GL;
typedef ProductQuadratureAP<AzimuthalU, PolarDGL>         Product_U_DGL;
typedef ProductQuadratureAP<AzimuthalU, PolarU>           Product_U_U;
typedef ProductQuadratureAP<AzimuthalU, PolarTY>          Product_U_TY;
//
typedef ProductQuadratureAP<AzimuthalS, PolarGL>          Product_S_GL;
typedef ProductQuadratureAP<AzimuthalS, PolarDGL>         Product_S_DGL;
typedef ProductQuadratureAP<AzimuthalS, PolarU>           Product_S_U;
typedef ProductQuadratureAP<AzimuthalS, PolarTY>          Product_S_TY;
//
typedef ProductQuadratureAP<AzimuthalASQR, PolarASDR>     Product_ASQR_ASDR;

} // end namespace detran_angle

#endif /* detran_angle_PRODUCTQUADRATUREAP_HH_ */

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.hh
//----------------------------------------------------------------------------//
