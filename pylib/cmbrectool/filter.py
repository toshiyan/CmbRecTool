
import numpy as np, healpy as hp
from analysis_class import *

ac = analysis_choice()
f  = filename()
r  = recfunc()

# sim
ac.sn = 10
cl = np.zeros((ac.sn,5,ac.oL[1]+1))

for i in range(ac.sn):
  salm = hp.fitsfunc.read_alm(f.salm[i+1])
  nalm = hp.fitsfunc.read_alm(f.nalm[i+1])
  cl[i,1,:] = hp.sphtfunc.alm2cl(salm,lmax=ac.oL[1])/r.w2
  cl[i,2,:] = hp.sphtfunc.alm2cl(nalm,lmax=ac.oL[1])/r.w2

mcl = np.average(cl,axis=0)
mcl[0,:] = r.L
mcl[3,:] = np.sqrt(mcl[1,:]*r.ilcl[1,:])
mcl[4,:] = mcl[3,:]/(mcl[1,:]+mcl[2,:])
np.savetxt(f.fl,mcl.T)

