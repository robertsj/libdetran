from detran_utilities import *
import numpy as np

print "inputdb: "
db = InputDB()
print dir(db)

db.put_int("number_groups", 2)
ng = db.get_int("number_groups")
assert(ng == 2)

db.put_dbl("keff", 1.123)
keff = db.get_dbl("keff")
assert(keff == 1.123)


v = vec_int(4, 0)
v[0] = 1
v[1] = 2
v[2] = 3
v[3] = 5
print v[3]
v2 = np.asarray(v)
print v2
print v
db.put_vec_int("tstvec", v)
v3 = db.get_vec_int("tstvec")

v4 = [0,1,2,3]
db.put_vec_int("tstvec2", v4)

vdbl = [0.123, 0.314, 0.785]
db.put_vec_dbl("dblvec", vdbl)
vdbl2 = db.get_vec_dbl("dblvec")


print "number_groups = ", ng
print "keff = ", keff
print np.asarray(v3)
print np.asarray(db.get_vec_int("tstvec2"))
print np.asarray(vdbl2)
