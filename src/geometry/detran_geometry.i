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

%include std_vec_typemap.i

%apply (std::vector<int> INPUTVECTOR, std::vector<int> INPUTVECTOR, std::vector<double> INPUTVECTOR, std::vector<double> INPUTVECTOR, std::vector<int> INPUTVECTOR)
       {(std::vector<int> xfm, std::vector<int> yfm, std::vector<double> xcme, std::vector<double> ycme, std::vector<int> mat_map)}

%apply (std::vector<int> INPUTVECTOR)
       {(vec_int xfm)}
%apply (std::vector<int> INPUTVECTOR)
       {(vec_int yfm)}
%apply (std::vector<int> INPUTVECTOR)
       {(vec_int mat_map)}

%apply (std::vector<double> INPUTVECTOR)
       {(vec_dbl xcme)}
%apply (std::vector<double> INPUTVECTOR)
       {(vec_dbl ycme)}

%rename(pyMesh2D) Mesh2D(vec_int xfm, vec_int yfm, vec_dbl xcme, vec_dbl ycme, vec_int mat_map);
/*%apply (std::vector<int> ARG_NAME)*/
/*       {(std::vector<int> xfm)}*/

/*%apply (std::vector<double> INPUTVECTOR)*/
/*       {(std::vector<double> value)}*/

%include "Definitions.hh"
//%include "Mesh.hh"
//%include "Mesh2D.hh"

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

namespace detran
{

class Mesh2D
{

public:

  Mesh2D(std::vector<int> xfm, std::vector<int> yfm, std::vector<double> xcme, std::vector<double> ycme, std::vector<int> mat_map);

};

} // end namespace detran_utils






