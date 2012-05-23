//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper3D.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper3D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER3D_HH_
#define SWEEPER3D_HH_

// Detran
#include "Sweeper.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Sweeper3D
 * \brief Sweeper for 3D problems.
 */
//---------------------------------------------------------------------------//
template <class EQ>
class Sweeper3D: public Sweeper<_3D>
{

public:

  typedef SP<Sweeper3D>                     SP_sweeper;
  typedef Sweeper<_3D>                      Base;
  //
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  //
  typedef Dimension<EQ::dimension>          D;
  typedef EQ                                Equation_T;
  typedef Boundary<_3D>                     Boundary_T;
  typedef typename Boundary_T::SP_boundary  SP_boundary;
  typedef typename
      BoundaryTraits<_3D>::value_type       boundary_flux_type;
  //
  typedef typename
      SweepSource<_3D>::SP_sweepsource      SP_sweepsource;
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
   */
  Sweeper3D(SP_input input,
            SP_mesh mesh,
            SP_material material,
            SP_quadrature quadrature,
            SP_state state,
            SP_boundary boundary,
            SP_sweepsource sweepsource)
  : Base(input,mesh,material,quadrature,
         state,boundary,sweepsource)
  , d_equation(mesh, material, quadrature, d_update_psi)
  {}

  /// Stop Eclipse warnings.
  virtual ~Sweeper3D(){}

  /*!
   *  \brief SP Constructor.
   */
  static detran::SP<Sweeper3D<EQ> >
  Create(detran::SP<InputDB>                    input,
         detran::SP<detran::Mesh>               mesh,
         detran::SP<detran::Material>           material,
         detran::SP<detran::Quadrature>         quadrature,
         detran::SP<detran::State>              state,
         detran::SP<detran::Boundary<_3D> >     boundary,
         detran::SP<detran::SweepSource<_3D> >  sweepsource)
  {
    SP_sweeper p;
    p = new Sweeper3D(input, mesh, material, quadrature,
                      state, boundary, sweepsource);
    return p;
  }

  /// Sweep.
  inline void sweep(moments_type &phi);

  /// Setup the equations for the group
  void setup_group(int g)
  {
    d_g = g;
    d_equation.setup_group(g);
  }

private:

  /// Equation
  EQ d_equation;

};

} // end namespace detran

#include "Sweeper3D.i.hh"

#endif /* SWEEPER3D_HH_ */
