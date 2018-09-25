
import numpy as np

x0 = .1 # mean and standard deviation
x = np.random.normal(0,x0,10000)
mc2 = np.average(x**2)
mc4 = np.average(x**4)
mc6 = np.average(x**6)
#sc4 = x**4 - mc**2
#hc4 = x**4 - 2*x**2*mc + mc**2
#print np.var(sc4), np.var(hc4)
#print np.mean(sc4), np.mean(hc4)

sc6 = -2*mc2**3 + 3*mc2**2*x**2
hc6 = (16./5.)*mc6
print np.var(sc6), np.var(hc6)
print np.mean(sc6), np.mean(hc6)
print mc6

