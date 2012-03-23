//----------------------------------*-C++-*----------------------------------//
/*!
 *  \file   Mesh3D.hh
 *  \author Jeremy Roberts
 *  \brief  Mesh3D class definition.
 */
//---------------------------------------------------------------------------//

// Geometry headers
#include "Mesh3D.hh"

namespace detran
{

Mesh3D::Mesh3D(vec_int xfm,  vec_int yfm,  vec_int zfm,
               vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
               vec_int mat_map)
  : Mesh(3,
         xfm,
         yfm,
         zfm,
         xcme,
         ycme,
         zcme,
         mat_map)
{
  /* ... */
}

Mesh3D::Mesh3D(vec_dbl xfme, vec_dbl yfme, vec_dbl zfme, vec_int mat_map)
  : Mesh(3,
         xfme,
         yfme,
         zfme,
         mat_map)
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Mesh3D.cc
//---------------------------------------------------------------------------//
