#!/usr/bin/env python
# Some code based on Sumner, et. a.l, JPhysD, 20 (1987), 1095-1101

from __future__ import division

u = 1e3
t = 0.001
R = [0.1, 0.11, 0.12]

n = len(R)
Sti = [u*t/(2*r) for r in R]

# Approximation - until figure out how to do nested sums
St = 1
for nn in xrange(n-1):
    St *= Sti[-1]*Sti[nn]*(1 - (R[nn]/R[nn+1])**2)
print St
