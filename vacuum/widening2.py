# same as the first one but with integrals
from sympy import *
from pylab import plot, show, title, legend
from scipy import zeros, linspace

s0 = 1
x0 = 2

D, l, x, s, b, k = symbols('Dlxsbk')
D = s + x * l
D = s + x * l

cc = k * D ** 3 / l

res = 1/simplify(integrate(1/cc/l,(l,0,b)))
inflim = limit(res, b, oo)

mn = 31
Llist = linspace(1,10,mn)
Con = zeros(mn)
Con2 = zeros(mn)
for m in range(0,mn):
    Con[m] = res.subs(b, Llist[m]).subs(x, x0).subs(s, s0).subs(k, 2e4*2.6e-4)

title('Integral')
plot(Llist, Con,label="Vacumm conductance")
if (x0 > 0):
    inflimn = inflim.subs(x, x0).subs(s, s0).subs(k, 2e4*2.6e-4)
    plot([Llist[0],Llist[-1]],[inflimn,inflimn],'--',label="Infinite length")
legend(loc='best')
show()
