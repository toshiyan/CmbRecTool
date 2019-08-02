
import numpy as np, healpy as hp

lmax = 3000
nside = 2048
Fl = np.ones(lmax)

calm, x = hp.sphtfunc.Alm.getlm(lmax)
print(calm)

for i in range(1):
  print(i)
  alm1 = hp.sphtfunc.almxfl(calm,Fl)
  alm2 = hp.sphtfunc.almxfl(calm,Fl)
  print('alm2map')
  map1 = hp.sphtfunc.alm2map(alm1,nside)
  map2 = hp.sphtfunc.alm2map(alm2,nside)
  cmap = map1*map2
  print('convolution')
  tlm  = hp.sphtfunc.map2alm(cmap,lmax=lmax)

