from detran import *

fm = [100]
cm = [0.0, 10.0]
mt = [0]
mat  = Material.Create(1,1,False)


# 1D
mesh1 = Mesh1D.Create(fm, cm, mt)
quad1 = GaussLegendre.Create(8)
acc1 = CMR(mesh1, mat, quad1)
acc1.initialize(2)
mesh_fine = acc1.get_mesh()
mesh_fine.display()
print mesh_fine.dimension()
mesh_coarse = acc1.get_coarse_mesh()
mesh_coarse.display()

print "fine to coarse:"
print acc1.fine_to_coarse(0, 0)
print acc1.fine_to_coarse(1, 0)
print acc1.fine_to_coarse(2, 0)

## 2D
#mesh2 = Mesh2D.Create(fm, fm, cm, cm, mt)
#quad2 = QuadrupleRange.Create(2)
#acc2 = CMR(mesh2, mat, quad2)
#acc2.initialize(7)
#mesh_fine = acc2.get_mesh()
#mesh_fine.display()
#print mesh_fine.dimension()
#mesh_coarse = acc2.get_coarse_mesh()
#mesh_coarse.display()


