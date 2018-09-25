
import numpy as np, healpy as hp
from analysis_class import *

ac = analysis_choice()
f  = filename()
r  = recfunc()

ac.sn = 10
Fl = np.loadtxt(f.fl,unpack=True)[4]
cl = np.zeros((ac.sn+1,ac.oL[1]+1))

#//// tau convolution ////#
for i in range(ac.sn+1):
  print i
  calm = hp.fitsfunc.read_alm(f.calm[i])
  alm1 = hp.sphtfunc.almxfl(calm,Fl*r.fil)
  alm2 = hp.sphtfunc.almxfl(calm,Fl*r.fil*r.ilcl[0,:])
  map1 = hp.sphtfunc.alm2map(alm1,ac.nside)
  map2 = hp.sphtfunc.alm2map(alm2,ac.nside)
  cmap = map1*map2
  tlm  = hp.sphtfunc.map2alm(cmap)
  hp.fitsfunc.write_alm(f.tlm[i],tlm)
  #cl[i,:] = hp.sphtfunc.alm2cl(tlm,lmax=ac.oL[1])/r.w4

#mcl = np.average(cl[1:,:],axis=0)
#mvl = np.var(cl[1:,:],axis=0)
#np.savetxt('cltt.dat',(np.array((r.L,cl[0,:],mcl,mvl))).T)

