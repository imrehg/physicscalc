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

xx = linspace(-3,3,5000)
dd = xx[1]-xx[0]
#yy = square(xx, 1.5)
#yy = sin(2*pi*xx*35) + sin(2*pi*xx*10+0.5)
yy = gauss(xx, 0.025)*sin(2*pi*xx*50)
dx= 0.2
#dx2 = 0.5
yy2 = 0.5*(gauss(xx+dx, 0.025) + gauss(xx-dx, 0.025))*sin(2*pi*xx*50)
#yy2 = 0.25*(gauss(xx+dx2, 0.025) + gauss(xx-dx2, 0.025) + gauss(xx+dx, 0.025) + gauss(xx-dx, 0.025))*sin(2*pi*xx*50)
#yy = gauss(xx, 0.2) + gauss(xx-1, 0.2) + gauss(xx+1, 0.2)

# ## Fourier Transform
X1 = fft(yy)
X2 = fft(yy2)
n = len(X1)
F1 = array(range(0,(n/2)))/n/dd

figure(1)
plot(F1, abs(X1[0:n/2]))
plot(F1, abs(X2[0:n/2]))

F = concatenate((F1, -1*F1[::-1]))
H = X2 / X1
# plot(F, real(H),'-')
# plot(F, imag(H),'.')


for nn in range(len(X1)):
    sf = 0.2
    X1[nn] *= cos(2*pi*sf*F[nn])


    # # Low pass filter
    # lim = 30
    # slope = 20
    # H[nn] = 1/(abs(F[nn]/lim)) if abs(F[nn])>lim else 1
    # X1[nn] *= H[nn]

#     # # gaussian Filter
#     # centre = 50
#     # lim = 3
#     # X1[nn] = X1[nn]*exp(-(abs(F[nn]-centre)/lim)**2)

#     # # # Sinc filter
#     # centre = 54
#     # lim = 1
#     # X1[nn] = X1[nn]*sin(pi*(abs(F[nn])-centre)/lim)/(pi*(abs(F[nn])-centre)/lim)

    # # Sharp bandpass filter
    # centre = 55
    # width = 0.5
    # # dp = 0.1
    # H[nn] = 0 if abs(abs(F[nn]) - centre) > width/2 else exp(-1j*pi/2)
    # X1[nn] *= H[nn] 

#     # Lorentz filter
#     # centre = 60
#     # g = 0.5
    
#     # X1[nn] = X1[nn]*g/((abs(F[nn])-centre)*2 + g**2/4)/pi/2

plot(F1, abs(X1[0:n/2]))
xlim(0,100)
xlabel('Frequency')
ylabel('|FFT|')
# #X = ifft(X1)

figure(2)
plot(xx, yy, label='Original')
plot(xx, yy2, label='Goal')
plot(xx, ifft(X1), label='Transformed')
#plot(xx, ifft(H), label='Transform function')
xlabel('time')
legend(loc='best')

show()
