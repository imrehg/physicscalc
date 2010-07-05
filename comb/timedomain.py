from __future__ import division
from numpy import *
from pylab import *

def square(x, T):
    out = array([])
    for t in x:
        w = 0 if abs(t) > T else 1
        out = append(out, w)
    return out

def gauss(x, s):
    return exp(-(x/s)**2/2)

xx = linspace(-1,1,5000)
dd = xx[1]-xx[0]
#yy = square(xx, 1.5)
#yy = sin(2*pi*xx*35) + sin(2*pi*xx*10+0.5)
yy = gauss(xx, 0.025)*sin(2*pi*xx*50)
#yy = gauss(xx, 0.2) + gauss(xx-1, 0.2) + gauss(xx+1, 0.2)

## Fourier Transform
X1 = fft(yy)
n = len(X1)
F1 = array(range(0,(n/2)))/n/dd

figure(1)
plot(F1, abs(X1[0:n/2]))

F = concatenate((F1, -1*F1[::-1]))

for nn in range(len(X1)):
    # # Low pass filter
    # lim = 55
    # slope = 100
    # X1[n] = 1/(slope*abs(F[n]/lim))*X1[n] if abs(F[n])>lim else X1[n]

    # # gaussian Filter
    # centre = 50
    # lim = 3
    # X1[nn] = X1[nn]*exp(-(abs(F[nn]-centre)/lim)**2)

    # # Sinc filter
    centre = 54
    lim = 1
    X1[nn] = X1[nn]*sin(pi*(abs(F[nn])-centre)/lim)/(pi*(abs(F[nn])-centre)/lim)

    # Lorentz filter
    # centre = 60
    # g = 0.5
    
    # X1[nn] = X1[nn]*g/((abs(F[nn])-centre)*2 + g**2/4)/pi/2

plot(F1, abs(X1[0:n/2]))
xlim(0,100)
xlabel('Frequency')
ylabel('|FFT|')
#X = ifft(X1)

figure(2)
plot(yy, label='Original')
plot(ifft(X1), label='Transformed')
xlabel('time')
legend(loc='best')
show()
