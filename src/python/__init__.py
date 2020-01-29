# src/python/detran/__init__.py
#
# Copyright (C) 2013 Jeremy Roberts <j.alyn.roberts@gmail.com>

from utilities import *
from angle import *
from callow import *
from orthog import *
from geometry import *
from material import *
from external_source import *
from boundary import *
from transport import *
from kinetics import *
from ioutils import *
from solvers import *
from postprocess import *

has_numpy = False
try :
  import numpy as np
  has_numpy = True
except ImportError :
  print "Warning: Could not import Numpy.  Numpy is highly recommended."
  print "         Detran utilities are unavailable without Numpy."
  
if has_numpy :
  try :
    from pydetranutils import *
  except ImportError :
    print "Warning: Could not import Detran utilities."

import sys

global __detran_python_manager__

# Instantiate this class to ensure it outlives Detran objects.  Then, its
# destructor can safely finalize PETsc and such.
class PyManager :
  def __init__(self, argv) :
    self.M = Manager()
    self.M.initialize(argv)

__detran_python_manager__ = PyManager(sys.argv)

import atexit
def pyfinalize():
  __detran_python_manager__.M.finalize()
atexit.register(pyfinalize)