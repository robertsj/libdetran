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
typedef ProductQuadratureAP<AzimuthalGL, PolarGC>         Product_GL_GC;
typedef ProductQuadratureAP<AzimuthalGL, PolarDGC>        Product_GL_DCL;
typedef ProductQuadratureAP<AzimuthalGL, PolarU>          Product_GL_U;
typedef ProductQuadratureAP<AzimuthalGL, PolarTY>         Product_GL_TY;
typedef ProductQuadratureAP<AzimuthalGL, PolarASDR>       Product_GL_ASDR;
//
typedef ProductQuadratureAP<AzimuthalDGL, PolarGL>        Product_DGL_GL;
typedef ProductQuadratureAP<AzimuthalDGL, PolarDGL>       Product_DGL_DGL;
typedef ProductQuadratureAP<AzimuthalDGL, PolarGC>        Product_DGL_GC;
typedef ProductQuadratureAP<AzimuthalDGL, PolarDGC>       Product_DGL_DCL;
typedef ProductQuadratureAP<AzimuthalDGL, PolarU>         Product_DGL_U;
typedef ProductQuadratureAP<AzimuthalDGL, PolarTY>        Product_DGL_TY;
typedef ProductQuadratureAP<AzimuthalDGL, PolarASDR>      Product_DGL_ASDR;
//
typedef ProductQuadratureAP<AzimuthalU, PolarGL>          Product_U_GL;
typedef ProductQuadratureAP<AzimuthalU, PolarDGL>         Product_U_DGL;
typedef ProductQuadratureAP<AzimuthalU, PolarGC>          Product_U_GC;
typedef ProductQuadratureAP<AzimuthalU, PolarDGC>         Product_U_DGC;
typedef ProductQuadratureAP<AzimuthalU, PolarU>           Product_U_U;
typedef ProductQuadratureAP<AzimuthalU, PolarTY>          Product_U_TY;
typedef ProductQuadratureAP<AzimuthalU, PolarASDR>        Product_U_ASDR;
typedef ProductQuadratureAP<AzimuthalU, PolarTG>          Product_U_TG;
//
typedef ProductQuadratureAP<AzimuthalS, PolarGL>          Product_S_GL;
typedef ProductQuadratureAP<AzimuthalS, PolarDGL>         Product_S_DGL;
typedef ProductQuadratureAP<AzimuthalS, PolarGC>          Product_S_GC;
typedef ProductQuadratureAP<AzimuthalS, PolarDGC>         Product_S_DGC;
typedef ProductQuadratureAP<AzimuthalS, PolarU>           Product_S_U;
typedef ProductQuadratureAP<AzimuthalS, PolarTY>          Product_S_TY;
typedef ProductQuadratureAP<AzimuthalS, PolarASDR>        Product_S_ASDR;
//
typedef ProductQuadratureAP<AzimuthalASQR, PolarGL>       Product_ASQR_GL;
typedef ProductQuadratureAP<AzimuthalASQR, PolarDGL>      Product_ASQR_DGL;
typedef ProductQuadratureAP<AzimuthalASQR, PolarGC>       Product_ASQR_GC;
typedef ProductQuadratureAP<AzimuthalASQR, PolarDGC>      Product_ASQR_DGC;
typedef ProductQuadratureAP<AzimuthalASQR, PolarU>        Product_ASQR_U;
typedef ProductQuadratureAP<AzimuthalASQR, PolarTY>       Product_ASQR_TY;
typedef ProductQuadratureAP<AzimuthalASQR, PolarASDR>     Product_ASQR_ASDR;
typedef ProductQuadratureAP<AzimuthalASQR, PolarTG>       Product_ASQR_TG;
} // end namespace detran_angle

#endif /* detran_angle_PRODUCTQUADRATUREAP_HH_ */

//----------------------------------------------------------------------------//
//              end of ProductQuadrature.hh
//----------------------------------------------------------------------------//
