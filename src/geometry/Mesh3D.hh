//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh3D.hh
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_MESH3D_HH_
#define detran_geometry_MESH3D_HH_

#include "Mesh.hh"
#ifdef DETRAN_ENABLE_BOOST
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#endif

namespace detran_geometry
{

//---------------------------------------------------------------------------//
/**
 *  @class Mesh3D
 *  @brief Three-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class GEOMETRY_EXPORT Mesh3D : public Mesh
{

public:

  typedef Mesh            Base;
  typedef Base::SP_mesh   SP_mesh;

  /**
   *  @brief Constructor.
   *
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    yfm         Fine meshes per coarse mesh in y dimension.
   *  @param    zfm         Fine meshes per coarse mesh in z dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    ycme        Coarse mesh edges y dimension.
   *  @param    zcme        Coarse mesh edges z dimension.
   *  @param    mat_map     Coarse mesh material map.
   */
  Mesh3D(vec_int xfm,  vec_int yfm,  vec_int zfm,
         vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
         vec_int mat_map);

  /**
   *  @brief Constructor.
   *
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    yfme        Fine mesh edges y dimension.
   *  @param    zfme        Fine mesh edges z dimension.
   *  @param    mat_map     Fine mesh material map.
   */
   Mesh3D(vec_dbl xfme, vec_dbl yfme, vec_dbl zfme, vec_int mat_map);

  /**
   *  @brief SP Constructor.
   *
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    yfm         Fine meshes per coarse mesh in y dimension.
   *  @param    zfm         Fine meshes per coarse mesh in z dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    ycme        Coarse mesh edges y dimension.
   *  @param    zcme        Coarse mesh edges z dimension.
   *  @param    mat_map     Coarse mesh material map.
   */
  static SP_mesh
  Create(vec_int xfm, vec_int yfm, vec_int zfm, vec_dbl xcme, vec_dbl ycme,
         vec_dbl zcme, vec_int mat_map)
  {
    SP_mesh p;
    p = new Mesh3D(xfm, yfm, zfm, xcme, ycme, zcme, mat_map);
    return p;
  }

  /**
   *  @brief SP Constructor.
   *
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    yfme        Fine mesh edges y dimension.
   *  @param    zfme        Fine mesh edges z dimension.
   *  @param    mat_map     Fine mesh material map.
   */
  static SP_mesh
  Create(vec_dbl xfme, vec_dbl yfme, vec_dbl zfme, vec_int mat_map)
  {
    SP_mesh p;
    p = new Mesh3D(xfme, yfme, zfme, mat_map);
    return p;
  }

protected:

  /**
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh3D() : Mesh(3) {}

private:

#ifdef DETRAN_ENABLE_BOOST
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
  }
#endif

};

} // end namespace detran_geometry

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(detran_geometry::Mesh3D)
#endif

#endif /* detran_geometry_MESH3D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh3D.hh
//---------------------------------------------------------------------------//
