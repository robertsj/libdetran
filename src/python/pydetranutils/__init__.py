# src/python/detran/pydetranutils/__init__.py
#
# Copyright (C) 2013 Jeremy Roberts <j.alyn.roberts@gmail.com>

try :
  print "importing detran mesh plotting utilities..."
  from mesh_plot import *
  print "importing detran quadrature plotting utilities..."
  from quad_plot import *
except :
  print("Error importing Detran plotting utilities")
