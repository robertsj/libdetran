//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mesh2D.hh
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MESH2D_HH_
#define MESH2D_HH_

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
class Mesh2D : public Mesh
{

public:

  typedef detran_utils::SP<Mesh2D>          SP_mesh;
  typedef Mesh                              Base;
  typedef Base::vec_int                     vec_int;
  typedef Base::vec_dbl                     vec_dbl;


  /*!
   *  \brief Constructor.
   *
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    yfm         Fine meshes per coarse mesh in y dimension.
   *  \param    xcme        Coarse mesh edges x dimension.
   *  \param    ycme        Coarse mesh edges y dimension.
   *  \param    mat_map     Coarse mesh material map.
   */
  Mesh2D(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    xfme        Fine mesh edges x dimension.
   *  \param    yfme        Fine mesh edges y dimension.
   *  \param    mat_map     Fine mesh material map.
   */
   Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int mat_map);

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfm         Fine meshes per coarse mesh in x dimension.
    *  \param    yfm         Fine meshes per coarse mesh in y dimension.
    *  \param    xcme        Coarse mesh edges x dimension.
    *  \param    ycme        Coarse mesh edges y dimension.
    *  \param    mat_map     Coarse mesh material map.
    */
   static detran_utils::SP<Mesh> Create(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map)
   {
     detran_utils::SP<Mesh> p;
     p = new Mesh2D(xfm, yfm, xcme, ycme, mat_map);
     return p;
   }

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfme        Fine mesh edges x dimension.
    *  \param    yfme        Fine mesh edges y dimension.
    *  \param    mat_map     Fine mesh material map.
    */
    static detran_utils::SP<Mesh> Create(vec_dbl xfme, vec_dbl yfme, vec_int mat_map)
    {
      detran_utils::SP<Mesh> p;
      p = new Mesh2D(xfme, yfme, mat_map);
      return p;
    }

protected:

  /*!
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh2D() : Mesh(2) {}

};

} // end namespace detran

#endif /* MESH2D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh2D.hh
//---------------------------------------------------------------------------//
