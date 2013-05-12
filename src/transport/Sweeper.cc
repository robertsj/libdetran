//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper.cc
 *  @author robertsj
 *  @date   Nov 7, 2012
 *  @brief  Sweeper member definitions.
 */
//---------------------------------------------------------------------------//

#include "Sweeper.t.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
Sweeper<D>::Sweeper(SP_input input,
                    SP_mesh mesh,
                    SP_material material,
                    SP_quadrature quadrature,
                    SP_state state,
                    SP_boundary boundary,
                    SP_sweepsource sweepsource)
  : d_input(input)
  , d_mesh(mesh)
  , d_material(material)
  , d_quadrature(quadrature)
  , d_state(state)
  , d_sweepsource(sweepsource)
  , d_update_psi(false)
  , d_adjoint(false)
  , d_number_sweeps(0)
  , d_update_boundary(false)
  , d_ordered_octants(std::pow((float)2, (int)D::dimension), 0)
{
  Require(d_input);
  Require(d_mesh);
  Require(d_material);
  Require(d_quadrature);
  Require(d_state);
  Require(boundary);
  Require(d_sweepsource);

  // Check whether we keep psi.
  if (d_input->check("store_angular_flux"))
    d_update_psi = (0 != d_input->get<int>("store_angular_flux"));

  // Check for the adjoint mode.
  if (d_input->check("adjoint"))
    d_adjoint = 0 != d_input->get<int>("adjoint");

  // Perform templated setup tasks.
  setup();

  // Setup the space-angle sweep indices
  setup_spatial_indices();
  setup_octant_indices(boundary);

}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::setup_group(const size_t g)
{
  d_g = g;
}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::set_update_psi(const bool v)
{
  d_update_psi = v;
}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::set_update_boundary(const bool v)
{
  d_update_boundary = v;
}

//---------------------------------------------------------------------------//
template <class D>
bool Sweeper<D>::update_boundary() const
{
  return d_update_boundary;
}

//---------------------------------------------------------------------------//
template <class D>
detran_utilities::size_t Sweeper<D>::number_sweeps() const
{
  return d_number_sweeps;
}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::set_adjoint(const bool adjoint)
{
  d_adjoint = adjoint;
  // Ensure the mesh sweeps are correct.
  setup_spatial_indices();
}

//---------------------------------------------------------------------------//
template <class D>
bool Sweeper<D>::is_adjoint() const
{
  return d_adjoint;
}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::set_tally(SP_tally tally)
{
  Require(tally);
  d_tally = tally;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::setup_spatial_indices()
{
  // Setup the spatial indexing ranges.  This eliminates an
  // indexing function.  Additionally, it automatically
  // chooses the right sweeping order for the adjoint problem.
  d_space_ranges.resize(d_quadrature->number_octants(),
                        vec2_int(D::dimension,
                                 vec_int(2, 0)));
  // Octants for each dimension that don't need sweep in negative direction
  int oct[3][4] = {{0,3,4,7}, {0,1,4,5}, {0,1,2,3}};
  for (size_t o = 0; o < d_quadrature->number_octants(); ++o)
  {
    for (size_t dim = 0; dim < D::dimension; ++dim)
    {
      if ((o == oct[dim][0] || o == oct[dim][1]  ||
           o == oct[dim][2] || o == oct[dim][3]   ))// && !d_adjoint)
      {
        d_space_ranges[o][dim][0] = 0;
        d_space_ranges[o][dim][1] = 1;
      }
      else
      {
        d_space_ranges[o][dim][0] = d_mesh->number_cells(dim) - 1;
        d_space_ranges[o][dim][1] = -1;
      }
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
void Sweeper<D>::setup_octant_indices(SP_boundary boundary)
{
  int no = d_quadrature->number_octants();

  // Initialize the ordered octants in a staggered pattern.  This
  // is the default for pure vacuum and pure reflection.
  for (int i = 0; i < no/2; ++i)
  {
    d_ordered_octants[2*i    ] = i;
    d_ordered_octants[2*i + 1] = no/2 + i;
  }

  // Sort the octants if their are mixed conditions
  if (boundary->has_reflective() && boundary->has_vacuum())
  {
    // Reset to cyclic for ordering.
    for (int i = 0; i < no; ++i)
      d_ordered_octants[i] = i;

    // Order the octants so that vacuum conditions start first
    vec_int count(no, 0);
    for (int side = 0; side < 2*D::dimension; side++)
      if (!boundary->is_reflective(side))
        for (int o = 0; o < no/2; o++)
          ++count[d_quadrature->incident_octant(side)[o]];

    // Sort the octants
    for (int i = 0; i < no; i++)
      for (int j = 0; j < no; j++)
        if (count[j] < count[i])
        {
          int o = d_ordered_octants[j];
          d_ordered_octants[j] = d_ordered_octants[i];
          d_ordered_octants[i] = o;
          o = count[j];
          count[j] = count[i];
          count[i] = o;
        }
  } // end reflective

  // Overwrite with user-specified order (useful to show effect of order)
  if (d_input->check("sweeper_ordered_octants"))
  {
    vec_int oo = d_input->get<vec_int>("sweeper_ordered_octants");
    Insist(oo.size() == no, "Wrong number of ordered octants.");
    d_ordered_octants = oo;
  }

  if (d_input->check("sweeper_print_octants"))
  {
    if (d_input->get<int>("sweeper_print_octants"))
    {
      std::cout << " SWEEPER ORDERED OCTANTS: ";
      for (int o = 0; o < no; o++) std::cout << " " << d_ordered_octants[o];
      std::cout << std::endl;
    }
  }


}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

TRANSPORT_INSTANTIATE_EXPORT(Sweeper<_1D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper<_2D>)
TRANSPORT_INSTANTIATE_EXPORT(Sweeper<_3D>)

} // end namespace detran
