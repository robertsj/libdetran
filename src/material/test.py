from detran_materials import *

print dir(Material)

# 2 groups, 1 material, 
mat = Material(2, 1, True) 

print "   number groups = ", mat.number_groups()
print "number materials = ", mat.number_materials()

mat.set_sigma_t(0, 0, 0.123)  # g0 total
mat.set_sigma_t(0, 1, 0.666)  # g1 total

mat.set_sigma_s(0, 0, 0, 1.0) # g0->g0
mat.set_sigma_s(0, 0, 1, 2.0) # g1->g0 etc.
mat.set_sigma_s(0, 1, 0, 3.0)
mat.set_sigma_s(0, 1, 1, 4.0)
mat.finalize()

print "scattering"
print mat.sigma_s(0, 0, 0)
print mat.sigma_s(0, 0, 1)
print mat.sigma_s(0, 1, 0)
print mat.sigma_s(0, 1, 1)
print "downscatter?"
print mat.downscatter()
print "pretty display?"
mat.display()