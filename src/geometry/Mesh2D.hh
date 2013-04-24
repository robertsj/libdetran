//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh2D.hh
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_MESH2D_HH_
#define detran_geometry_MESH2D_HH_

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
 *  @class Mesh2D
 *  @brief Two-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class GEOMETRY_EXPORT Mesh2D : public Mesh
{

public:

  typedef Mesh            Base;
  typedef Base::SP_mesh   SP_mesh;

  /**
   *  @brief Constructor.
   *
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    yfm         Fine meshes per coarse mesh in y dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    ycme        Coarse mesh edges y dimension.
   *  @param    mat_map     Coarse mesh material map.
   */
  Mesh2D(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map);

  /**
   *  @brief Constructor.
   *
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    yfme        Fine mesh edges y dimension.
   *  @param    mat_map     Fine mesh material map.
   */
  Mesh2D(vec_dbl xfme, vec_dbl yfme, vec_int mat_map);

  /// SP_constructor
  static SP_mesh
  Create(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map)
  {
    SP_mesh p(new Mesh2D(xfm, yfm, xcme, ycme, mat_map));
    return p;
  }

  /// SP constructor
  static SP_mesh
  Create(vec_dbl xfme, vec_dbl yfme, vec_int mat_map)
  {
    SP_mesh p(new Mesh2D(xfme, yfme, mat_map));
    return p;
  }

protected:

  /**
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh2D() : Mesh(2) {}

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

} // end namespace detran

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(detran_geometry::Mesh2D)
#endif

#endif /* detran_geometry_MESH2D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh2D.hh
//---------------------------------------------------------------------------//
