//----------------------------------*-C++-*----------------------------------//
/*!
 *  \file   Mesh.cc
 *  \author Jeremy Roberts
 *  \brief  Mesh class member definitions.
 */
//---------------------------------------------------------------------------//

// Geometry headers
#include "Mesh.hh"

// System headers
#include <numeric>
#include <iostream>

namespace detran
{

Mesh::Mesh(int dim,
           vec_int xfm,  vec_int yfm,  vec_int zfm,
           vec_dbl xcme, vec_dbl ycme, vec_dbl zcme,
           vec_int mat_map)
  : d_dimension(dim)
  , d_xfm(xfm)
  , d_yfm(yfm)
  , d_zfm(zfm)
  , d_xcme(xcme)
  , d_ycme(ycme)
  , d_zcme(zcme)
{
  Require(mat_map.size() == d_xfm.size()*d_yfm.size()*d_zfm.size());
  // Setup discretizations, etc.
  setup();
  // Set the coarse mesh material map.
  std::string s = "MATERIAL";
  add_coarse_mesh_map(s, mat_map);
};

Mesh::Mesh(int dim,
           vec_dbl xfme, vec_dbl yfme, vec_dbl zfme,
           vec_int mat_map)
  : d_dimension(dim)
  , d_xfm(vec_int(xfme.size()-1, 1))
  , d_yfm(vec_int(yfme.size()-1, 1))
  , d_zfm(vec_int(zfme.size()-1, 1))
  , d_xcme(xfme)
  , d_ycme(yfme)
  , d_zcme(zfme)
{
  // Setup discretizations, etc.
  setup();
  // Add the fine mesh material map.
  Assert(mat_map.size() == d_number_cells);
  std::string s = "MATERIAL";
  add_mesh_map(s, mat_map);
};



void Mesh::add_coarse_mesh_map(std::string map_key, vec_int m_map)
{
  // Temporary map.
  vec_int tmp_m_map(d_number_cells, 0);

  int i_save = 0; // place holders
  int j_save = 0;
  int k_save = 0;

  for (int k = 0; k < d_zfm.size(); k++)
  {
    // Fine mesh z range for this kth coarse mesh.
    int k1 = k_save;
    int k2 = k_save + d_zfm[k];

    for (int j = 0; j < d_yfm.size(); j++)
    {
      // Fine mesh x range for this jth coarse mesh.
      int j1 = j_save;
      int j2 = j_save + d_yfm[j];

      for (int i = 0; i < d_xfm.size(); i++)
      {
        // Fine mesh y range for this ith coarse mesh.
        int i1 = i_save;
        int i2 = i_save + d_xfm[i];

        for (int kk = k1; kk < k2; kk++)
        {
          for (int jj = j1; jj < j2; jj++)
          {
            for (int ii = i1; ii < i2; ii++)
            {
              tmp_m_map[ii + jj * d_number_cells_x  +
                        kk * d_number_cells_x * d_number_cells_y] =
                  m_map[i + j * d_xfm.size() + k * d_xfm.size()*d_yfm.size()];
            }
          }
        }
        i_save = std::accumulate(d_xfm.begin(), d_xfm.begin() + i + 1, 0);
      }
      i_save = 0;
      j_save = std::accumulate(d_xfm.begin(), d_xfm.begin() + j + 1, 0);
    }
    j_save = 0;
    k_save = std::accumulate(d_zfm.begin(), d_zfm.begin() + k + 1, 0);
  }

  add_mesh_map(map_key, tmp_m_map);
  return;
}

/*!
 *
 */
void Mesh::add_mesh_map(std::string map_key, vec_int mesh_map)
{
  Require(!map_key.empty());
  Require(mesh_map.size() == d_number_cells);

  // Erase the value associated with the key if it exists.
  mesh_map_type::iterator it;
  it = d_mesh_map.find(map_key);
  if (it != d_mesh_map.end())
    d_mesh_map.erase(it);

  // Add the new value.
  d_mesh_map[map_key] = mesh_map;

}

// Check if fine mesh map exists.
bool Mesh::mesh_map_exists(std::string map_key)
{
  mesh_map_type::iterator iter;
  iter = d_mesh_map.find(map_key);
  if (iter != d_mesh_map.end())
    return true;
  else
    return false;
}

/*!
 *
 */
const vec_int& Mesh::mesh_map(std::string map_key)
{
  // Add the new value.
  return d_mesh_map[map_key];
}


void Mesh::setup()
{
  Require(d_xfm.size() > 0);
  Require(d_yfm.size() > 0);
  Require(d_zfm.size() > 0);
  Require(d_xcme.size() == d_xfm.size()+1);
  Require(d_ycme.size() == d_yfm.size()+1);
  Require(d_zcme.size() == d_zfm.size()+1);

  // Compute numbers of cells.
  d_number_cells_x = std::accumulate(d_xfm.begin(), d_xfm.end(), 0);
  d_number_cells_y = std::accumulate(d_yfm.begin(), d_yfm.end(), 0);
  d_number_cells_z = std::accumulate(d_zfm.begin(), d_zfm.end(), 0);
  d_number_cells   =  d_number_cells_x * d_number_cells_y * d_number_cells_z;

  // Cell widths.
  d_dx.resize(d_number_cells_x, 0.0);
  d_dy.resize(d_number_cells_y, 0.0);
  d_dz.resize(d_number_cells_z, 0.0);

  // Total domain widths
  d_total_width_x = d_xcme[d_xcme.size()-1] - d_xcme[0];
  d_total_width_y = d_ycme[d_ycme.size()-1] - d_ycme[0];
  d_total_width_z = d_ycme[d_zcme.size()-1] - d_zcme[0];

  // Discretize.
  int ph = 0; // place holder
  for (int i = 0; i < d_xfm.size(); i++)
  {
    // Fine mesh x range for this Ith coarse mesh.
    int i1 = ph;
    int i2 = ph + d_xfm[i];
    for (int ii = i1; ii < i2; ii++)
      d_dx[ii] = (d_xcme[i + 1] - d_xcme[i]) / d_xfm[i];
    ph = std::accumulate(d_xfm.begin(), d_xfm.begin() + i + 1, 0);

  }
  ph = 0;
  for (int j = 0; j < d_yfm.size(); j++)
  {
    // Fine mesh y range for this Jth coarse mesh.
    int j1 = ph;
    int j2 = ph + d_yfm[j];
    for (int jj = j1; jj < j2; jj++)
    {
      Assert(jj < d_number_cells_y);
      d_dy[jj] = (d_ycme[j + 1] - d_ycme[j]) / d_yfm[j];
    }
    ph = std::accumulate(d_yfm.begin(), d_yfm.begin() + j + 1, 0);
  }
  ph = 0;
  for (int k = 0; k < d_zfm.size(); k++)
  {
    // Fine mesh z range for this Kth coarse mesh.
    int k1 = ph;
    int k2 = ph + d_zfm[k];
    for (int kk = k1; kk < k2; kk++)
      d_dz[kk] = (d_zcme[k + 1] - d_zcme[k]) / d_zfm[k];
    ph = std::accumulate(d_zfm.begin(), d_zfm.begin() + k + 1, 0);
  }

}

void Mesh::display() const
{
  using std::cout;
  using std::endl;
  cout << endl << "Detran Mesh" << endl;
  cout << "           dimension: " << d_dimension << endl;
  cout << "        number cells: " << d_number_cells << endl;
  cout << "      number x cells: " << d_number_cells_x << endl;
  cout << "      number y cells: " << d_number_cells_y << endl;
  cout << "      number z cells: " << d_number_cells_z << endl;
  cout << "             x width: " << d_total_width_x << endl;
  cout << "             y width: " << d_total_width_y << endl;
  cout << "             z width: " << d_total_width_z << endl;
  cout << " x coarse mesh edges: " << endl << "   ";
  for (int i = 0; i < d_xcme.size(); i++)
  {
    cout << d_xcme[i] << " ";
  }
  cout << endl << " y coarse mesh edges: " << endl << "   ";
  for (int i = 0; i < d_ycme.size(); i++)
  {
    cout << d_ycme[i] << " ";
  }
  cout << endl << " z coarse mesh edges: " << endl << "   ";
  for (int i = 0; i < d_zcme.size(); i++)
  {
    cout << d_zcme[i] << " ";
  }
  cout << endl << " x fine mesh count: " << endl << "   ";
  for (int i = 0; i < d_xfm.size(); i++)
  {
    cout << d_xfm[i] << " ";
  }
  cout << endl << " y fine mesh count: " << endl << "   ";
  for (int i = 0; i < d_yfm.size(); i++)
  {
    cout << d_yfm[i] << " ";
  }
  cout << endl << " z fine mesh count: " << endl << "   ";
  for (int i = 0; i < d_zfm.size(); i++)
  {
    cout << d_zfm[i] << " ";
  }
  cout << endl << endl;
}



} // end namespace detran

//---------------------------------------------------------------------------//
//              end of Mesh.cc
//---------------------------------------------------------------------------//
