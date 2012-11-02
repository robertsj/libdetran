//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.hh
 * \author Jeremy Roberts
 * @date   Mar 24, 2012
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
  typedef typename Base::SP_tally                   SP_tally;
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
   *  @param    input       User input database.
   *  @param    mesh        Cartesian mesh.
   *  @param    material    Material database.
   *  @param    quadrature  Angular quadrature.
   *  @param    state       State vectors.
   *  @param    boundary    Boundary based on mesh.
   *  @param    sweepsource Sweep source constructor.
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
  , d_ordered_octants(4, 0)
  {
    // Default order - staggered.  I've notices that
    // doing cyclic sweeps (0->1->2->3) leads to anisotropic
    // edge fluxes when they should be isotropic (thought the
    // cell centered psi is fine)
    d_ordered_octants[0] = 0;
    d_ordered_octants[1] = 2;
    d_ordered_octants[2] = 1;
    d_ordered_octants[3] = 3;

    // Order the octants so that vacuum conditions start first
    vec_int count(4, 0);
    for (int side = 0; side < 4; side++)
    {
      if (!d_boundary->is_reflective(side))
      {
        for (int o = 0; o < 2; o++)
          ++count[d_quadrature->incident_octant(side)[o]];
      }
    }
    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        if (count[j] < count[i])
        {
          int o = d_ordered_octants[j];
          d_ordered_octants[j] = d_ordered_octants[i];
          d_ordered_octants[i] = o;
          o = count[j];
          count[j] = count[i];
          count[i] = o;
        }
      }
    }
    std::cout << " ORDERED OCTANTS: " << std::endl;
    for (int o = 0; o < 4; o++)
      std::cout << " o = " << d_ordered_octants[o] << std::endl;

  }

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
  vec_int d_ordered_octants;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "Sweeper2D.i.hh"

#endif /* SWEEPER2D_HH_ */
