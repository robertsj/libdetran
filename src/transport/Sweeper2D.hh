//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER2D_HH_
#define SWEEPER2D_HH_

#include "Sweeper.hh"
#include "boundary/BoundarySN.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Sweeper2D
 *  \brief Sweeper for 2D discrete ordinates problems.
 */
//---------------------------------------------------------------------------//
template <class EQ>
class Sweeper2D: public Sweeper<_2D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Sweeper2D>           SP_sweeper;
  typedef Sweeper<_2D>                              Base;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::Mesh                       Mesh;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_quadrature              SP_quadrature;
  typedef typename Base::SP_sweepsource             SP_sweepsource;
  typedef typename Base::moments_type               moments_type;
  typedef typename Base::angular_flux_type          angular_flux_type;
  typedef typename Base::SP_currenttally            SP_currenttally;
  typedef typename Base::vec_int                    vec_int;
  typedef typename Base::vec2_int                   vec2_int;
  typedef typename Base::vec3_int                   vec3_int;
  typedef typename Base::size_t                     size_t;
  typedef EQ                                        Equation_T;
  typedef BoundarySN<_2D>                           Boundary_T;
  typedef typename Boundary_T::SP_boundary          SP_boundary;
  typedef typename BoundaryTraits<_2D>::value_type  boundary_flux_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   *  \param    boundary    Boundary based on mesh.
   *  \param    sweepsource Sweep source constructor.
   */
  Sweeper2D(SP_input input,
            SP_mesh mesh,
            SP_material material,
            SP_quadrature quadrature,
            SP_state state,
            SP_boundary boundary,
            SP_sweepsource sweepsource)
  : Base(input,mesh,material,quadrature,
         state,sweepsource)
  , d_boundary(boundary)
  {}

  /// Virtual destructor
  virtual ~Sweeper2D(){}

  /// SP Constructor
  static SP_sweeper
  Create(SP_input       input,
         SP_mesh        mesh,
         SP_material    material,
         SP_quadrature  quadrature,
         SP_state       state,
         SP_boundary    boundary,
         SP_sweepsource sweepsource)
  {
    SP_sweeper p(new Sweeper2D(input, mesh, material, quadrature,
                               state, boundary, sweepsource));
    return p;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SWEEPERS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Sweep.
  inline void sweep(moments_type &phi);

  /// Setup the equations for the group
  void setup_group(const size_t g)
  {
    d_g = g;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_boundary d_boundary;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "Sweeper2D.i.hh"

#endif /* SWEEPER2D_HH_ */
