//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ProductQuadrature.t.cc
 *  @brief ProductQuadrature explicit instantiations
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//----------------------------------------------------------------------------//

#include "angle/ProductQuadratureAP.t.hh"

namespace detran_angle
{

// product
//
template class ProductQuadratureAP<AzimuthalGL, PolarGL>;
template class ProductQuadratureAP<AzimuthalGL, PolarDGL>;
template class ProductQuadratureAP<AzimuthalGL, PolarU>;
template class ProductQuadratureAP<AzimuthalGL, PolarTY>;
//
template class ProductQuadratureAP<AzimuthalDGL, PolarGL>;
template class ProductQuadratureAP<AzimuthalDGL, PolarDGL>;
template class ProductQuadratureAP<AzimuthalDGL, PolarU>;
template class ProductQuadratureAP<AzimuthalDGL, PolarTY>;
//
template class ProductQuadratureAP<AzimuthalU, PolarGL>;
template class ProductQuadratureAP<AzimuthalU, PolarDGL>;
template class ProductQuadratureAP<AzimuthalU, PolarU>;
template class ProductQuadratureAP<AzimuthalU, PolarTY>;
//
template class ProductQuadratureAP<AzimuthalS, PolarGL>;
template class ProductQuadratureAP<AzimuthalS, PolarDGL>;
template class ProductQuadratureAP<AzimuthalS, PolarU>;
template class ProductQuadratureAP<AzimuthalS, PolarTY>;
//
template class ProductQuadratureAP<AzimuthalASQR, PolarASDR>;


} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of ProductQuadratureAP.t.cc
//----------------------------------------------------------------------------//
