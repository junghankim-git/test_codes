#!/usr/bin/env python

import matplotlib.pyplot as plt


total = 15475.995
initial = 3.364+166.995
dynamics = 11453.027
physics = 1609.861
write = 48.752+1076.042
others = total-initial-dynamics-physics-write


print 'total = ', total
print 'initial = ', initial
print 'dynamics = ', dynamics
print 'physics = ', physics
print 'write  = ', write
print 'others = ', others
print 'total = ', 100.
print 'initial = ', initial/total*100.
print 'dynamics = ', dynamics/total*100.
print 'physics = ', physics/total*100.
print 'write  = ', write/total*100.
print 'others = ', others/total*100.
