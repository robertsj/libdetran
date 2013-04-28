# src/python/detran/pydetranutils/__init__.py
#
# Copyright (C) 2013 Jeremy Roberts <j.alyn.roberts@gmail.com>

try :
  from mesh_plot import *
  from quad_plot import *
except :
  print("Error importing Detran plotting utilities")
