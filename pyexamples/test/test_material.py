from detran import Material

# Simple 1 group, 1 material

mat = Material(1, 1, True)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_s(0, 0, 0, 0.5)
mat.finalize()
assert( mat.sigma_t(0, 0)    == 1.0 )
assert( mat.sigma_s(0, 0, 0) == 0.5 )

# Two group with upscatter

mat = Material(2, 1, False)
mat.set_sigma_t(0, 0,    1.0)
mat.set_sigma_t(0, 1,    2.0)
mat.set_sigma_s(0, 0, 0, 0.3)
mat.set_sigma_s(0, 0, 1, 0.1)
mat.set_sigma_s(0, 1, 0, 0.2)
mat.set_sigma_s(0, 1, 1, 0.9)
mat.finalize()
mat.display()
assert( mat.sigma_t(0, 0)    == 1.0 )
assert( mat.sigma_s(0, 1, 0) == 0.2 )
