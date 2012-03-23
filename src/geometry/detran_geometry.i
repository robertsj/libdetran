//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_geometry.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

%module detran_geometry
%{
#include "Definitions.hh"
#include "Mesh.hh"
#include "Mesh2D.hh"
%}

// Load the standard library interfaces
%include std_vector.i
%include std_string.i

// Load the vector maps.  Note, I am a bit unhappy with
// how it's all used.  They work *if* I declare the class
// interface below.  Otherwise, just including e.g.
// Mesh2D doesn't allow the maps, since I'm using 
// typedefs on the input arguments.  There should be an
// easy way around this, but I'm not a SWIG pro.
%include std_vec_typemap.i

%apply (std::vector<int>    INPUTVECTOR, 
        std::vector<int>    INPUTVECTOR, 
        std::vector<double> INPUTVECTOR, 
        std::vector<double> INPUTVECTOR, 
        std::vector<int>    INPUTVECTOR)
      {(std::vector<int>    xfm, 
        std::vector<int>    yfm, 
        std::vector<double> xcme, 
        std::vector<double> ycme, 
        std::vector<int>    mat_map)}

%include "Definitions.hh"
%include "Mesh.hh"

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

namespace detran
{
// Simple redefinition of the basic interface.
class Mesh2D : public Mesh
{
public:
  Mesh2D(std::vector<int>    xfm, 
         std::vector<int>    yfm, 
         std::vector<double> xcme, 
         std::vector<double> ycme, 
         std::vector<int>    mat_map);
  Mesh2D(std::vector<double> xfme, 
         std::vector<double> yfme, 
         std::vector<int>    mat_map);
public:
  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
  int index(int i, int j = 0, int k = 0);
  //virtual void mesh_map(std::string map_key, std::vector<int> mesh_map);
};


} // end namespace detran_utils





