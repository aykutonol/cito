#!/usr/bin/env python3

from __future__ import print_function

# Mass
m = 1.
# Half size
# red box: 0.15075 0.0745 0.0398
lx = 0.15075
ly = 0.0745
lz = 0.0398
# Full size
Lx = lx*2
Ly = ly*2
Lz = lz*2
# Calculate diagonal inertia terms
Ixx = 1/12*m*(Ly**2 + Lz**2)
Iyy = 1/12*m*(Lx**2 + Lz**2)
Izz = 1/12*m*(Lx**2 + Ly**2)
# Print the result
print("Mass: ", m)
print("Lx: %f, Ly: %f, Lz: %f" % (Lx, Ly, Lz))
print("Ixx: %.12f, Iyy: %.12f, Izz: %.12f" % (Ixx, Iyy, Izz))
print("Copy and paste to diaginertia:\n%.12f %.12f %.12f" % 
      (Ixx, Iyy, Izz))