//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh3D.hh
 *  @author Jeremy Roberts
 *  @brief  Mesh3D member definitions.
 */
//---------------------------------------------------------------------------//

#include "Mesh3D.hh"

namespace detran_geometry
{

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
Mesh3D::Mesh3D(vec_dbl xfme, vec_dbl yfme, vec_dbl zfme, vec_int mat_map)
  : Mesh(3,
         xfme,
         yfme,
         zfme,
         mat_map)
{
  /* ... */
}

} // end namespace detran_geometry

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(detran_geometry::Mesh3D)
#endif

//---------------------------------------------------------------------------//
//              end of Mesh3D.cc
//---------------------------------------------------------------------------//
