from detran_materials import *

# 2 groups, 1 material, 
mat = Material(2, 1, True) 

print "   number groups = ", mat.number_groups()
print "number materials = ", mat.number_materials()

mat.set_sigma_t(0, 0, 1.0)    # g0 total
mat.set_sigma_t(0, 1, 1.0)    # g1 total

mat.set_sigma_s(0, 0, 0, 1.0) # g0->g0
mat.set_sigma_s(0, 0, 1, 1.0) # g1->g0 etc.
mat.set_sigma_s(0, 1, 0, 1.0)
mat.set_sigma_s(0, 1, 1, 1.0)

mat.finalize()

print mat.downscatter()