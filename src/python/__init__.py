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

from pydetranutils import *

import numpy as np
import sys
import atexit

print "*** detran -- python interface ***"

Manager.initialize(sys.argv)

def goodbye() :
  print "*** te ootte nyt loppuun ***"
  Manager.finalize()
atexit.register(goodbye)

