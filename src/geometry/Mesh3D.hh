//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mesh3D.hh
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MESH3D_HH_
#define MESH3D_HH_

// Geometry headers
#include "Mesh.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Mesh2D
 *  \brief Two-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class Mesh3D : public Mesh
{

public:

  typedef detran_utils::SP<Mesh3D>          SP_mesh;
  typedef Mesh                              Base;
  typedef Base::vec_int                     vec_int;
  typedef Base::vec_dbl                     vec_dbl;


  /*!
   *  \brief Constructor.
   *
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    yfm         Fine meshes per coarse mesh in y dimension.
   *  \param    zfm         Fine meshes per coarse mesh in z dimension.
   *  \param    xcme        Coarse mesh edges x dimension.
   *  \param    ycme        Coarse mesh edges y dimension.
   *  \param    zcme        Coarse mesh edges z dimension.
   *  \param    mat_map     Coarse mesh material map.
   */
  Mesh3D(vec_int xfm,  vec_int yfm,  vec_int zfm,
         vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
         vec_int mat_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    xfme        Fine mesh edges x dimension.
   *  \param    yfme        Fine mesh edges y dimension.
   *  \param    zfme        Fine mesh edges z dimension.
   *  \param    mat_map     Fine mesh material map.
   */
   Mesh3D(vec_dbl xfme, vec_dbl yfme, vec_dbl zfme, vec_int mat_map);

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfm         Fine meshes per coarse mesh in x dimension.
    *  \param    yfm         Fine meshes per coarse mesh in y dimension.
    *  \param    zfm         Fine meshes per coarse mesh in z dimension.
    *  \param    xcme        Coarse mesh edges x dimension.
    *  \param    ycme        Coarse mesh edges y dimension.
    *  \param    zcme        Coarse mesh edges z dimension.
    *  \param    mat_map     Coarse mesh material map.
    */
   static SP_mesh CreateMesh3D(vec_int xfm,  vec_int yfm,  vec_int zfm,
                               vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
                               vec_int mat_map)
   {
     SP_mesh p;
     p = new Mesh3D(xfm, yfm, zfm, xcme, ycme, zcme, mat_map);
     return p;
   }

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfme        Fine mesh edges x dimension.
    *  \param    yfme        Fine mesh edges y dimension.
    *  \param    zfme        Fine mesh edges z dimension.
    *  \param    mat_map     Fine mesh material map.
    */
    static SP_mesh CreateMesh3D(vec_dbl xfme,
                                vec_dbl yfme,
                                vec_dbl zfme,
                                vec_int mat_map)
    {
      SP_mesh p;
      p = new Mesh3D(xfme, yfme, zfme, mat_map);
      return p;
    }


protected:

  /*!
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh3D() : Mesh(3) {}

};

} // end namespace detran


#endif /* MESH3D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh3D.hh
//---------------------------------------------------------------------------//
