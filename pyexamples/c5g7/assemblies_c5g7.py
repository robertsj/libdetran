# pyexamples/assembly.py
#
#,2-D pin cell
import numpy as np
import time
from detran import *

def get_pins() :

  # Shared
  pitch = 1.26
  radii = [0.54]
  number = 7
  flag = True
  # Pin 0 - UO2 
  pin0 = PinCell.Create(pitch, radii, [0,6])
  pin0.meshify(number, flag)
  # Pin 1 - 4.3% MOX
  pin1 = PinCell.Create(pitch, radii, [1,6])
  pin1.meshify(number, flag)
  # Pin 2 - 7.0% MOX
  pin2 = PinCell.Create(pitch, radii, [2,6])
  pin2.meshify(number, flag)
  # Pin 3 - 8.7% MOX
  pin3 = PinCell.Create(pitch, radii, [3,6])
  pin3.meshify(number, flag)
  # Pin 4 - Guide Tube
  pin4 = PinCell.Create(pitch, radii, [4,6])
  pin4.meshify(number, flag)
  # Pin 5 - Fission Chamber
  pin5 = PinCell.Create(pitch, radii, [5,6])
  pin5.meshify(number, flag)
  # Pin 6 - Moderator
  pin6 = PinCell.Create(pitch, radii, [6,6])
  pin6.meshify(number, flag)

  return pin0, pin1, pin2, pin3, pin4, pin5, pin6

def get_assemblies() :

  # Shared things
  G = 4 # guide tube
  F = 5 # fission chamber
  pin0, pin1, pin2, pin3, pin4, pin5, pin6 = get_pins()

  # Assembly 1 -- UO2
  assem1 = Assembly.Create(0)
  assem1.add_pincell(pin0)
  assem1.add_pincell(pin1)
  assem1.add_pincell(pin2)
  assem1.add_pincell(pin3)
  assem1.add_pincell(pin4)
  assem1.add_pincell(pin5)
  assem1.add_pincell(pin6)
  pin_map1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
              0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,F,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,G,0,0,G,0,0,G,0,0,G,0,0,G,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,G,0,0,0,0,0,0,0,0,0,G,0,0,0,
              0,0,0,0,0,G,0,0,G,0,0,G,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  assem1.finalize(pin_map1)

  # Assembly 2 - MOX
  assem2 = Assembly.Create(17)
  assem2.add_pincell(pin0)
  assem2.add_pincell(pin1)
  assem2.add_pincell(pin2)
  assem2.add_pincell(pin3)
  assem2.add_pincell(pin4)
  assem2.add_pincell(pin5)
  assem2.add_pincell(pin6)
  pin_map2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
              1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
              1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
              1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
              1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
              1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,G,3,3,G,3,3,F,3,3,G,3,3,G,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,1,
              1,2,G,3,3,G,3,3,G,3,3,G,3,3,G,2,1,
              1,2,2,2,3,3,3,3,3,3,3,3,3,2,2,2,1,
              1,2,2,G,2,3,3,3,3,3,3,3,2,G,2,2,1,
              1,2,2,2,2,G,2,2,G,2,2,G,2,2,2,2,1,
              1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,
              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
  assem2.finalize(pin_map2)

  # Assembly 3 - Moderator
  assem3 = Assembly.Create(17)
  assem3.add_pincell(pin0)
  assem3.add_pincell(pin1)
  assem3.add_pincell(pin2)
  assem3.add_pincell(pin3)
  assem3.add_pincell(pin4)
  assem3.add_pincell(pin5)
  assem3.add_pincell(pin6)
  pin_map3 = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
              6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];
  assem3.finalize(pin_map3)

  return assem1, assem2, assem3


