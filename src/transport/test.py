from detran_transport import *
from detran_angle import *
from detran_geometry import *
from detran_materials import *
from detran_utilities import *

#print "transport test"
#print dir(InputDB)

# Make input
input = InputDB.Create()
input.put_int("number_groups", 2)
input.put_int("store_angular_flux", 1)

# Make mesh
xcm = [0.0, 1.0, 2.0]
xfm = [  10, 10   ]
ycm = [0.0, 1.0, 2.0]
yfm = [  10, 10   ]
mat_map = [ 0,  0, 
            0,  0 ]
mesh = Mesh2D.Create(xfm, yfm, xcm, ycm, mat_map)
#print mesh

# Make quadrature
quadrature = QuadrupleRange.Create(2)
print quadrature
print dir(quadrature)

quadrature.display()

# Make state
print input
print input.this
print mesh
print mesh.this
state = State.Create(input.this, mesh.this, quadrature.this)
print state

