//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputDB.c
 * \author Jeremy Roberts
 * \brief  InputDB class definition.
 */
//---------------------------------------------------------------------------//


#include "InputDB.hh"

#include <iostream>

namespace detran
{

InputDB::InputDB()
{
  /* ... */
}

void InputDB::display() const
{
  using std::cout;
  using std::endl;
  using std::string;
  cout << "InputDB contents" << endl;
  {
    // Integers
    cout << "  Integer data: " << endl;
    std::map<string, int>::const_iterator it = d_data_int.begin();
    for (; it != d_data_int.end(); it++)
    {
      cout << "    " << (*it).first << "  " << (*it).second << endl;
    }
  }
  {
    // Doubles
    cout << "  Double data: " << endl;
    std::map<string, double>::const_iterator it = d_data_dbl.begin();
    for (; it != d_data_dbl.end(); it++)
    {
      cout << "    " << (*it).first << "  " << (*it).second << endl;
    }
  }
  {
    // Strings
    cout << "  String data: " << endl;
    std::map<string, string>::const_iterator it = d_data_str.begin();
    for (; it != d_data_str.end(); it++)
    {
      cout << "    " << (*it).first << "  " << (*it).second << endl;
    }
  }
  {
    // Integer vectors
    cout << "  Integer vector data: " << endl;
    std::map<string, vec_int>::const_iterator it = d_data_vec_int.begin();
    for (; it != d_data_vec_int.end(); it++)
    {
      cout << "    " << (*it).first << "  ";
      vec_int v = (*it).second;
      for (int i = 0; i < v.size(); i++) cout << v[i] << " ";
      cout << endl;
    }
  }
  {
    // Double vectors
    cout << "  Double vector data: " << endl;
    std::map<string, vec_dbl>::const_iterator it = d_data_vec_dbl.begin();
    for (; it != d_data_vec_dbl.end(); it++)
    {
      cout << "    " << (*it).first << "  ";
      vec_dbl v = (*it).second;
      for (int i = 0; i < v.size(); i++) cout << v[i] << " ";
      cout << endl;
    }
  }
}

int InputDB::size(int type) const
{
  int value = 0;
  if (type == INT)
    value =  d_data_int.size();
  else if (type == DBL)
    value = d_data_dbl.size();
  else if (type == STR)
    value = d_data_str.size();
  else if (type == VEC_INT)
    value = d_data_vec_int.size();
  else if (type == VEC_DBL)
    value = d_data_vec_dbl.size();
  else
    THROW("Invalid type specifier.");
  return value;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of InputDB.cc
//---------------------------------------------------------------------------//
