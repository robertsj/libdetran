# /pyexamples/test/ioutils

from detran import *

def test_IO_HDF5() :
  """ Test the IO_HDF5 Python interface.
  """

  # Open the hdf5 file
  io = IO_HDF5("test.h5")
  
  # Create the test input to be written.
  inp = InputDB.Create()
  inp.put_int("number_groups", 2)
  inp.put_dbl("initial_keff",  1.234)
  inp.put_str("bc_left",       "vacuum")
  inp.put_vec_int("testvi",    [0,1,2,3])
  inp.put_vec_dbl("testvd",    [2.1, 3.4])

  # Create the test material to be written.
  mat = Material.Create(2, 2, False)
  mat.set_sigma_t(0, 0, 0.1890);             
  mat.set_sigma_s(0, 0, 0, 0.1507);       
  mat.set_sigma_f(1, 0, 0.0067); 
  mat.set_sigma_s(1, 1, 1, 0.9355);  

  # Write out and close.
  io.write(inp)
  io.write(mat)
  io.close()

  # Fill a new input 
  inp2 = io.read_input()
  assert(inp2.get_int("number_groups") == 2)
  assert(inp2.get_dbl("initial_keff")  == 1.234)
  
  # Fill a new material 
  mat2 = io.read_material()
  assert(mat2.sigma_t(0, 0) == 0.1890)

  return True

if __name__ == "__main__" :

  if test_IO_HDF5() :
    print "test_IO_HDF5 passed"

