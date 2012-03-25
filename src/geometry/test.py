from detran_geometry import *
import numpy as np
print "geometry test"
print dir(Mesh2D)

# Simple test.  A 2x2 coarse grid with 3x3 fine
# meshes per coarse grid.

# Fine mesh resolution.
N = 3 

# First make define via python lists.
xcm = [0.0, 1.0, 2.0]
xfm = [    N,   N   ]
ycm = [0.0, 1.0, 2.0]
yfm = [    N,   N   ]
# A material map and some other test map.  Note, the C++
# routines DO NOT flip the maps.  Moreover, the maps
# are all row-oriented.  That means the following map
# as formatted represents the map physically left-to-right 
# and top-to-bottom.  Python helper functions could be
# added to manage this, but that's not needed for now.
m_map_py = [ 0,  2,   
             1,  3 ]  
t_map_py = [ 1,  3,   
             6,  9]   

# Then via C++ vectors
cm = vec_dbl(3, 0.0)
fm = vec_int(2, N)
cm[0] = 0.0
cm[1] = 1.0
cm[2] = 2.0
m_map_cc    = vec_int(4, 0)
m_map_cc[0] = 0
m_map_cc[1] = 2
m_map_cc[2] = 1
m_map_cc[3] = 3
t_map_cc    = vec_int(4, 0)
t_map_cc[0] = 1
t_map_cc[1] = 3
t_map_cc[2] = 6
t_map_cc[3] = 9

# with python
mesh_py = Mesh2D(xfm, yfm, xcm, ycm, m_map_py)
mesh_py_sp = Mesh2D.CreateMesh2D(xfm, yfm, xcm, ycm, m_map_py)

# with c++
mesh_cc = Mesh2D(fm, fm, cm, cm, m_map_cc)

# short constructor.
#mesh_py2 = Mesh2D(cm, cm, m_map_py)

print "are these the same? : "
print mesh_py.dx(0),            mesh_cc.dx(0), mesh_py_sp.dx(0)
print mesh_py.number_cells(),   mesh_cc.number_cells()
print mesh_py.number_cells_x(), mesh_cc.number_cells_x()
print mesh_py.number_cells_y(), mesh_cc.number_cells_y()

print "dx"
for i in range(0, mesh_py.number_cells_x()) :
    print mesh_py.dx(i)
    assert(mesh_py.dx(i) == 1.0/N)

# Add maps
mesh_py.add_coarse_mesh_map("TMPMAP", t_map_py)
mesh_cc.add_coarse_mesh_map("TMPMAP", t_map_cc)

t_map_py2 = mesh_py.mesh_map("TMPMAP")
t_map_cc2 = mesh_cc.mesh_map("TMPMAP")
print np.asarray(t_map_py2)
print np.asarray(t_map_cc2)

# Plots
mesh_py.plot_mesh_map("MATERIAL")
mesh_py.plot_mesh_map("TMPMAP")