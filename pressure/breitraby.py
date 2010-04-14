
from __future__ import division
from scipy import sign, array, power
from pylab import plot, show

def breitrabi(params):
    Ahfs, I, gi, gj, B, mi, mj = params
    # ub = 9.27400899e-24
    ub =  1.399624624e6
    dEhfs = Ahfs * (I + 1/2)
    m = mi + mj
    s = sign(mj)
    x = (gj - gi) * ub * B / dEhfs
    p1 = - dEhfs / (2 * (2 * I + 1))
    p2 = gi * ub * m * B
    p3 = s * dEhfs / 2 * (1 + (4 * m * x)/(2 * I + 1) + power(x,2)) ** (1/2)
    if (m == -I-abs(mj)) & (mi < 0):
        b = 4 * m / (2 * I + 1)
        z = - b / 2
        for i in range(len(p3)):
            if (x[i] < z):
                p3[i] *= -1
    return p1 + p2 + p3
#    return s * (1 + (4 * m * x) / (2 * I + 1) + x**2)


h = 6.62606876e-34
Ahfs = 2.2981579425e9
I = 7 / 2
gi = -0.00039885395
gj = 2.00254032
B = array(range(0,10000,1))
#mi = 7/2
#mj = 1/2
#params = (Ahfs, I, gi, gj, B, mi, mj)
#plot(B, breitrabi(params)/1e9)

J = 1/2
for F in [3, 4]:
    for mf in range(-F, F+1):
        mj = F - I
        mi = mf - mj
        if (abs(mi) > I):
            mi = sign(mi)*(I)
            mj *= -1
        if mj >= 0:
            r = 'k'
        else:
            r = 'r'
        params = (Ahfs, I, gi, gj, B, mi, mj)
        plot(B, breitrabi(params)/1e9, r)

show()
