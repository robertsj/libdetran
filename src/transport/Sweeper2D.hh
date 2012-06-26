//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER2D_HH_
#define SWEEPER2D_HH_

// Detran
#include "Sweeper.hh"
#include "Boundary.hh"

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

  typedef SP<Sweeper2D>                     SP_sweeper;
  typedef Sweeper<_2D>                      Base;
  //
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  //
  typedef EQ                                Equation_T;
  typedef Boundary<_2D>                     Boundary_T;
  typedef typename Boundary_T::SP_boundary  SP_boundary;
  typedef typename
      BoundaryTraits<_2D>::value_type       boundary_flux_type;
  //
  typedef typename
      SweepSource<_2D>::SP_sweepsource      SP_sweepsource;
  //
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

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
  static detran::SP<Sweeper2D<EQ> >
  Create(detran::SP<InputDB>                    input,
         detran::SP<detran::Mesh>               mesh,
         detran::SP<detran::Material>           material,
         detran::SP<detran::Quadrature>         quadrature,
         detran::SP<detran::State>              state,
         detran::SP<detran::BoundaryBase<_2D> > boundary,
         detran::SP<detran::SweepSource<_2D> >  sweepsource)
  {
    SP_sweeper p(new Sweeper2D(input, mesh, material, quadrature,
                               state, boundary, sweepsource));
    return p;
  }

  /// Sweep.
  inline void sweep(moments_type &phi);

  /// Setup the equations for the group
  void setup_group(int g)
  {
    d_g = g;
  }

private:

  SP_boundary d_boundary;

};

} // end namespace detran

#include "Sweeper2D.i.hh"

#endif /* SWEEPER2D_HH_ */
