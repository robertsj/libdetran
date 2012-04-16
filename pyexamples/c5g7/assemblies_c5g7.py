# pyexamples/assembly.py
#
# 2-D pin cell
import numpy as np
import time
from detran import *

def get_pins() :

  # Shared
  pitch = 1.26
  radii = [0.54]
  number = 20
  # Pin 1 - UO2 
  pin1 = PinCell.Create(pitch, radii, [0,6])
  pin1.meshify(number)
  # Pin 2 - 4.3% MOX
  pin2 = PinCell.Create(pitch, radii, [1,6])
  pin2.meshify(number)
  # Pin 3 - 7.0% MOX
  pin3 = PinCell.Create(pitch, radii, [2,6])
  pin3.meshify(number)
  # Pin 4 - 8.7% MOX
  pin4 = PinCell.Create(pitch, radii, [3,6])
  pin4.meshify(number)
  # Pin 5 - Guide Tube
  pin5 = PinCell.Create(pitch, radii, [4,6])
  pin5.meshify(number)
  # Pin 6 - Fission Chamber
  pin6 = PinCell.Create(pitch, radii, [5,6])
  pin6.meshify(number)
  # Pin 7 - Moderator
  pin7 = PinCell.Create(pitch, [], [6])
  pin7.meshify(number)

  return pin1, pin2, pin3, pin4, pin5, pin6, pin7

def get_assemblies() :
  # Assembly 1 -- UO2
  assem1 = Assembly.Create(17)
  pin1, pin2, pin3, pin4, pin5, pin6, pin7 = get_pins()
  assem1.add_pincell(pin1)
  assem1.add_pincell(pin2)
  assem1.add_pincell(pin3)
  assem1.add_pincell(pin4)
  assem1.add_pincell(pin5)
  assem1.add_pincell(pin6)
  assem1.add_pincell(pin7)
  G = 4
  F = 5
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
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];     
  assem1.finalize(pin_map1)
  return assem1


