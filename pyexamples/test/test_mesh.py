from detran import Mesh1D, Mesh2D, Mesh3D

cm = [0.0, 2.0, 4.0]
fm = [   10,  10   ]

# 1D Mesh
mt1 = [0, 1]
mesh1 = Mesh1D.Create(fm, cm, mt1)
mesh1.display()

# 2D Mesh
mt2 = [0, 1, \
       0, 1]
mesh2 = Mesh2D.Create(fm, fm, cm, cm, mt2)
mesh2.display()

# 3D Mesh
mt3 = [0, 1, \
       0, 1, \
         1, 0, \
         1, 0]
mesh3 = Mesh3D.Create(fm, fm, fm, cm, cm, cm, mt3)
mesh3.display()
