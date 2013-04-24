//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearExternalSource.i.hh
 *  @brief  LinearExternalSource.i
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_LINEAREXTERNALSOURCE_I_HH_
#define detran_LINEAREXTERNALSOURCE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline double LinearExternalSource::source(const size_t cell,
                                           const size_t group)
{
  // Preconditions
  Require(cell  < d_mesh->number_cells());
  Require(group < d_number_groups);

  double value = 0.0;
  value += d_fa * d_sources[d_ia]->source(cell, group);
  value += d_fb * d_sources[d_ib]->source(cell, group);
  return value;
}

//---------------------------------------------------------------------------//
inline double LinearExternalSource::source(const size_t cell,
                                           const size_t group,
                                           const size_t angle)
{
  // Preconditions
  Require(cell  < d_mesh->number_cells());
  Require(group < d_number_groups);
  Require(angle < d_number_angles);

  double value = 0.0;
  value += d_fa * d_sources[d_ia]->source(cell, group, angle);
  value += d_fb * d_sources[d_ib]->source(cell, group, angle);
  return value;
}

//---------------------------------------------------------------------------//
inline void LinearExternalSource::set_time(const double time)
{
  d_time = time;

  // Determine the indices and factors, i.e.
  //   source = fa*sources[ia] + fb*sources[ib]
  d_ia = 0;
  d_ib = 0;
  d_fa = 1.0;
  d_fb = 0.0;
  if (d_time <= d_times[0])
  {
    /* ... */
  }
  else if (d_time > d_times[d_number_times - 1])
  {
    d_ia = d_number_times - 1;
    d_ib = d_ia;
  }
  else
  {
    for (int i = 1; i < d_number_times; ++i)
      if (d_time > d_times[i - 1] && d_time <= d_times[i])
        d_ib = i;
    d_ia = d_ib - 1;
    d_fb = (d_time - d_times[d_ia]) / (d_times[d_ib] - d_times[d_ia]);
    d_fa = 1.0 - d_fb;
  }

}

} // end namespace detran

#endif // LINEAREXTERNALSOURCE_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file LinearExternalSource.i.hh
//---------------------------------------------------------------------------//
