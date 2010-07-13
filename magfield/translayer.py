#!/usr/bin/env python
# Some code based on Sumner, et. a.l, JPhysD, 20 (1987), 1095-1101

from __future__ import division

def partsum(start, end, sums, level, S, R):
    """ Recursive function for partial sum
    Gives the different terms of Equation (11) in the paper
    start,end: summation limits
    sums: how many summation terms are there
    level: which sum is being calculated (0: outermost, ...)
    S, R: individual transverse screening factors and radii
    """
    if (level == sums):
        return 1
    finish = end - sums + level + 1
    out = 0
    for nn in range(start, finish):
        if level == 0:
            out += S[nn]*partsum(nn+1, end, sums, level+1, S, R)
        else:
            out += S[nn]*(1-(R[start-1]/R[nn])**2)*partsum(nn+1, end, sums, level+1, S, R)
    return out

def getscreening(Sti, R):
    """ Do the whole Equation (11)
    Sti, R: individual screening factors and radii
    """
    n = len(Sti)
    out = 1
    for sums in range(1, n+1):
        temp = partsum(0, n, sums, 0,  Sti, R)
        out += temp
    return out

if __name__ == "__main__":
    u = 1e3
    t = 0.001
    R = [0.1, 0.12, 0.17]
    Sti = [u*t/(2*r) for r in R]
    n = len(Sti)
    print "Transverse screening factor: %f" %(getscreening(Sti, R))
