# src/python/detran/pydetranutils/__init__.py
#
# Copyright (C) 2013 Jeremy Roberts <j.alyn.roberts@gmail.com>

global __detranuselatex__ 
__detranuselatex__ = False

try :
  import mpl_toolkits.mplot3d.axes3d as p3
  import matplotlib.pyplot as plt
except ImportError :
  print "Warning: Could not import matplotlib."

from mesh_plot import *
from quad_plot import *