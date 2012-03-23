from detran_geometry import *

print "geometry test"
print dir(Mesh2D)
xcm = [0.0, 1.0, 2.0]
xfm = [  10, 10   ]
ycm = [0.0, 1.0, 2.0]
yfm = [  10, 10   ]
mat_map = [ 0,  0, 
            0,  0 ]


cm = vec_dbl(3, 0.0)
fm = vec_int(2, 10)
cm[0] = 0.0
cm[1] = 3.0
cm[2] = 6.0
mm = vec_int(4, 0)

mesh = Mesh2D(fm, fm, cm, cm, mm)

#print mesh.dx(0)
#print mesh.number_cells()
#print mesh.number_cells_x()
#print mesh.dimension()
mesh2 = Mesh2D(fm, yfm, xcm, ycm, mat_map)
print mesh2.dx(0)

mesh3 = Mesh2D(xcm, ycm, mat_map)
print mesh3.dx(0)