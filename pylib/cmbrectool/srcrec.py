
import numpy as np, healpy as hp, iof90
from scipy.io import FortranFile

def nconv(alm,Lmax,ns):
  map = hp.sphtfunc.alm2map(alm,ns)
  return hp.sphtfunc.map2alm(map**2,lmax=lmax)

dir = '../../data/'
w = 'Lmask1'
lmin  = 100
lmax  = 2048
nside = 2048
simn  = 1

cl = (np.loadtxt(dir+'aps/aps_sim_1d_smica_'+w+'.dat',usecols=(1,2,3))).T
bl = (np.loadtxt(dir+'public/beam_smica.dat')).T
fl = np.zeros(lmax+1)
gl = np.zeros(lmax+1)
L  = np.linspace(0,lmax,lmax+1)
fl[lmin:] = 1./(cl[0,lmin-1:lmax]+cl[1,lmin-1:lmax]+cl[2,lmin-1:lmax])/bl[lmin-1:lmax]
gl[lmin:] = 1./(cl[0,lmin-1:lmax]+cl[1,lmin-1:lmax]+cl[2,lmin-1:lmax])

mask = iof90.readf90map(dir+'map/mask_'+w+'.bin')

cl = np.zeros((simn+1,lmax+1))
omap = iof90.readf90map(dir+'map/real_smica.bin')
oalm = hp.sphtfunc.almxfl(hp.sphtfunc.map2alm(omap*mask,lmax=lmax),fl)
alm = nconv(oalm,lmax,nside)
cl[0,:] = hp.sphtfunc.alm2cl(alm,lmax=lmax)

for i in range(1,simn+1):
  smap = iof90.readf90map(dir+'map/lcdm_smica_'+str(i)+'.bin')
  nmap = iof90.readf90map(dir+'map/nois_smica_'+str(i)+'.bin')
  pmap = iof90.readf90map(dir+'map/ptsr_smica_'+w+'_'+str(i)+'.bin')
  salm = hp.sphtfunc.almxfl(hp.sphtfunc.map2alm(smap*mask,lmax=lmax),fl)
  nalm = hp.sphtfunc.almxfl(hp.sphtfunc.map2alm(nmap*mask,lmax=lmax),fl)
  palm = hp.sphtfunc.almxfl(hp.sphtfunc.map2alm(pmap*mask,lmax=lmax),gl)
  alm = nconv(salm+nalm+palm,lmax,nside)
  cl[i,:] = hp.sphtfunc.alm2cl(alm,lmax=lmax)

mcl = np.average(cl[1:,:],axis=0)
mvl = np.sqrt(np.var(cl[1:,:],axis=0))
np.savetxt('clnn.dat',(np.array((L,cl[0,:],mcl,mvl))).T)

